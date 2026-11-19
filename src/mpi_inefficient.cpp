#include <mpi.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <algebra/finite_fields/finite_fields.hpp>
#include <algebra/matrices/matrices.hpp>

namespace {

using Field = algebra::FiniteField<11, 3>;
using Element = Field::Element;

constexpr std::size_t dimension = ALGEBRA_ASSIGNMENT_SIZE;
constexpr std::uint32_t seed = 123456u;
using Matrix = algebra::Matrix<Element, dimension, dimension>;

static_assert(algebra::FieldElement<Element>);
static_assert(algebra::RingElement<Matrix>);

double seconds(double start) {
    return MPI_Wtime() - start;
}

std::uint32_t draw_coefficient(std::mt19937& generator, std::uint32_t prime) {
    if (prime < 2) {
        throw std::invalid_argument("The prime must be at least two");
    }
    constexpr std::uint64_t range = std::uint64_t{1} << 32;
    const std::uint64_t limit = range - range % prime;
    for (;;) {
        const std::uint64_t value = generator();
        if (value < limit) {
            return static_cast<std::uint32_t>(value % prime);
        }
    }
}

Element next_element(std::mt19937& generator) {
    Element::polynomial_type polynomial;
    for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        const auto coefficient = draw_coefficient(generator, Field::characteristic);
        polynomial.set_coefficient(degree, Element::coefficient_type::from_integer(coefficient));
    }
    return Element::from_polynomial(polynomial);
}

void populate(Matrix& a, Matrix& b) {
    std::mt19937 generator(seed);
    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            a(i, j) = next_element(generator);
            b(i, j) = next_element(generator);
        }
    }
}

Element determinant(Matrix matrix, MPI_Comm comm) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Element result = Element::one();
    for (std::size_t k = 0; k < dimension; ++k) {
	// don't care about the pivot decision
        std::size_t pivot_row = k;
        while (pivot_row < dimension && matrix(pivot_row, k).is_zero()) {
            ++pivot_row;
        } // not even caring about anything but a non-zero pivot
        if (pivot_row == dimension) {
            return Element::zero();
        }
        if (pivot_row != k) {
            matrix.swap_rows(pivot_row, k); // actual swap
            result = result.negate(); // negation as we take product
        }
        const Element pivot = matrix(k, k); // the element by which we are pivoting
        result = result.multiply(pivot); // very singular operation
        Element inverse = Element::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        // this is the low-hanging fruit to parallelize
	// each worker should only calculate its own rows
	std::vector<uint64_t> buffer; // create a LOCAL contiguous buffer we can broadcast across MPI
        for (std::size_t i = k + 1 + static_cast<size_t>(rank); i < dimension; i += static_cast<size_t>(size)) {
            if (matrix(i, k).is_zero()) {
                for (size_t j = k + 1; j < dimension; ++j) { // pack our buffer no matter what
			buffer.push_back(Element::get_index(matrix(i, j).representative()));
		}
                continue;
            }
            const Element multiplier = matrix(i, k).multiply(inverse);
            matrix(i, k) = Element::zero();

            for (std::size_t j = k + 1; j < dimension; ++j) {
                const Element correction = multiplier.multiply(matrix(k, j));
                matrix(i, j) = matrix(i, j).add(correction.negate());
	    	buffer.push_back(Element::get_index(matrix(i, j).representative())); // further maintain contiguous buffer in single row
            }
        }

	// create a GLOBAL contiguous buffer that can be receieved via MPI
	std::vector<int> recv_n(size);
	std::vector<int> recv_disp(size); // create a log of displacements to tell MPI where to place data in global buffer
	int count = 0;
	for (int r = 0; r < size; r++) { // iterate each rank
		size_t owned = 0;
		for (size_t i = k + 1 + static_cast<size_t>(r); i < dimension; i += static_cast<size_t>(size)) owned++;
		recv_n[r] = static_cast<int>(owned * (dimension - (k + 1))); // true size of LOCAL contiguous buffer in GLOBAL contiguous buffer for each rank
		recv_disp[r] = count; // where THIS rank's data belongs in GLOBAL contiguous buffer
		count += recv_n[r]; // total size of GLOBAL contiguous buffer (all data to be broadcast)
	}
	std::vector<uint64_t> recv_buffer(count); // allocate GLOBAL contiguous buffer
	MPI_Allgatherv(
		buffer.data(),
		static_cast<int>(buffer.size()),
		MPI_UINT64_T,
		recv_buffer.data(),
		recv_n.data(),
		recv_disp.data(),
		MPI_UINT64_T,
		comm
	);
	for (int r = 0; r < size; r++) { // iterate each rank
		if (r == rank) continue; // the LOCAL rank already did this by hand
		size_t loc = static_cast<size_t>(recv_disp[r]); // access the real location of current rank's data in GLOBAL buffer
		for (size_t i = k + 1 + static_cast<size_t>(r); i < dimension; i += static_cast<size_t>(size)) {
			matrix(i, k) = Element::zero();
			for (size_t j = k + 1; j < dimension; j++) {
				matrix(i, j) = Element::from_index(recv_buffer[loc++]); // copy it over
			}
		}
	}
	// everything synchronized after MPI distribution so move to the next pivot
    }
    return result;
}

void print_determinant(const char* name, const Element& value, double elapsed) {
    std::cout << name << " = " << value << '\n';
    std::cout << "Coefficients [";
    for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        if (degree != 0) {
            std::cout << ", ";
        }
        std::cout << value.coefficient(degree);
    }
    std::cout << "]\nSeconds " << elapsed << '\n' << std::flush;
}

}

int main(int argc, char** argv) {
    try {
	MPI_Init(&argc, &argv);
	int rank = 0;
	int size = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

        const auto total_start = MPI_Wtime();
        Field::initialize();

	if (rank == 0)
        std::cout << std::fixed << std::setprecision(6)
                  << "Matrix dimension " << dimension << '\n'
                  << "Seed " << seed << '\n'
                  << "Field order " << Field::cardinality() << '\n'
                  << "Modulus " << Field::modulus() << '\n';

	if (rank == 0)
        std::cout << "Generating A and B\n" << std::flush;
        auto start = MPI_Wtime();
        Matrix a;
        Matrix b;
        populate(a, b);
	if (rank == 0)
        std::cout << "Generation seconds " << seconds(start) << '\n';

	if (rank == 0)
        std::cout << "Computing C = A * B\n" << std::flush;
        start = MPI_Wtime();
        const Matrix c = a.multiply(b);
	if (rank == 0)
        std::cout << "Product seconds " << seconds(start) << '\n';

	if (rank == 0)
        std::cout << "Computing det(A)\n" << std::flush;
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const auto det_a = determinant(a, MPI_COMM_WORLD);
	if (rank == 0)
        print_determinant("det(A)", det_a, seconds(start));

	if (rank == 0)
        std::cout << "Computing det(B)\n" << std::flush;
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const auto det_b = determinant(b, MPI_COMM_WORLD);
	if (rank == 0)
        print_determinant("det(B)", det_b, seconds(start));

	if (rank == 0)
        std::cout << "Computing det(C)\n" << std::flush;
        MPI_Barrier(MPI_COMM_WORLD);s
        start = MPI_Wtime();
        const auto det_c = determinant(c, MPI_COMM_WORLD);
	if (rank == 0)
        print_determinant("det(C)", det_c, seconds(start));

        const auto expected = det_a.multiply(det_b);
        const bool verified = det_c == expected;
	if (rank == 0)
        std::cout << "det(A) * det(B) = " << expected << '\n'
                  << "Identity verified " << std::boolalpha << verified << '\n'
                  << "Total seconds " << seconds(total_start) << '\n';

	MPI_Finalize();
        return verified ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
