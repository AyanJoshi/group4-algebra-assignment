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
	/* // WE NEED TO TRANSMIT LESS DATA
	// don't care about the pivot decision
        std::size_t pivot_row = k;
        while (pivot_row < dimension && matrix(pivot_row, k).is_zero()) {
            ++pivot_row; // don't care about anything but a non-zero pivot
        }
	*/
	uint64_t pivot_idx = std::numeric_limits<uint64_t>::max(); // default "nothing" state
	for (size_t candidate = static_cast<size_t>(rank); candidate < dimension; candidate += static_cast<size_t>(size)) {
		if (candidate < k) continue;
		if (!matrix(candidate, k).is_zero()) {
			pivot_idx = static_cast<uint64_t>(candidate); // take the first non-zero pivot
			break;
		}
	}
	uint64_t recv_pivot = std::numeric_limits<uint64_t>::max(); // default "nothing" state
	
	// make sure the pivot is the same on every rank
	MPI_Allreduce(
		&pivot_idx,
		&recv_pivot,
		1,
		MPI_UINT64_T,
		MPI_MIN,
		comm
	);

	/*
        if (pivot_row == dimension) {
            return Element::zero();
        }
	*/
	if (recv_pivot == std::numeric_limits<uint64_t>::max()) return Element::zero(); // same semantic singularity check

	/*
        if (pivot_row != k) {
            matrix.swap_rows(pivot_row, k); // actual swap
            result = result.negate(); // negation as we take product
        }
	*/
	std::vector<uint64_t> row; // create a LOCAL contiguous buffer we can broadcast across MPI
	int recv_owner = static_cast<int>(recv_pivot % static_cast<size_t>(size));
	int k_owner = static_cast<int>(k % static_cast<size_t>(size));
        if (recv_pivot != k) { // if we are actually doing anything
		if (k_owner == recv_owner && rank == k_owner) { // then just swap
			matrix.swap_rows(recv_pivot, k);
		} else if (rank == recv_owner) { // then send the pivot we're swapping!
			for (size_t i = k; i < dimension; i++) { // populate row buffer from owner rank
				row.push_back(Element::get_index(matrix(recv_pivot, i).representative()));
			}
			MPI_Sendrecv_replace(
				row.data(),
				static_cast<int>(row.size()),
				MPI_UINT64_T,
				k_owner,
				0,
				k_owner,
				0,
				comm,
				MPI_STATUS_IGNORE
			);
			matrix.swap_rows(recv_pivot, k); // go ahead and swap anyway so we don't overwrite our real data
			for (size_t i = 0; i < dimension - k; i++) {
				matrix(recv_pivot, i + k) = Element::from_index(row[i]);
			}
		} else if (rank == k_owner) { // then send the k we're swapping!
			for (size_t i = k; i < dimension; i++) { // populate row buffer from owner rank
				row.push_back(Element::get_index(matrix(k, i).representative()));
			}
			MPI_Sendrecv_replace(
				row.data(),
				static_cast<int>(row.size()),
				MPI_UINT64_T,
				recv_owner,
				0,
				recv_owner,
				0,
				comm,
				MPI_STATUS_IGNORE
			);
			matrix.swap_rows(recv_pivot, k); // go ahead and swap anyway so we don't overwrite our real data
			for (size_t i = 0; i < dimension - k; i++) {
				matrix(k, i + k) = Element::from_index(row[i]);
			}
		}
        	result = result.negate(); // negation as we take product
	}
	
	// now that pivot is authoritative, we need the pivot populated on ranks that had no authority over the swap
	std::vector<uint64_t> final_row(dimension - k); // create a LOCAL contiguous buffer we can broadcast across MPI
	if (rank == k_owner) { // take from a rank with a real pivot
		for (size_t i = k; i < dimension; i++) {
			final_row[i - k] = Element::get_index(matrix(k, i).representative()); // pivot data is now in row k
		}
	}

	MPI_Bcast( // it would be best if we didn't redo some of the send/recv work up above but it doesn't matter that much
		final_row.data(),
		static_cast<int>(final_row.size()),
		MPI_UINT64_T,
		k_owner,
		comm
	);
	
	if (rank != k_owner) { // overwrite the pivot to k on rows that got the copy (again, duplicates some work)
		for (size_t i = 0; i < dimension - k; i++) {
			matrix(k, i + k) = Element::from_index(final_row[i]);
		}
	}

        const Element pivot = matrix(k, k); // the element by which we are pivoting
        result = result.multiply(pivot); // very singular operation
        Element inverse = Element::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        for (std::size_t i = static_cast<size_t>(rank); i < dimension; i += static_cast<size_t>(size)) {
	    if (i <= k) continue;
            if (matrix(i, k).is_zero()) continue;
            const Element multiplier = matrix(i, k).multiply(inverse);
            matrix(i, k) = Element::zero();

            for (std::size_t j = k + 1; j < dimension; ++j) {
                const Element correction = multiplier.multiply(matrix(k, j));
                matrix(i, j) = matrix(i, j).add(correction.negate());
            }
        }
	// we don't care about real synchronization
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
        start = MPI_Wtime();
        const auto det_a = determinant(a, MPI_COMM_WORLD);
	if (rank == 0)
        print_determinant("det(A)", det_a, seconds(start));

	if (rank == 0)
        std::cout << "Computing det(B)\n" << std::flush;
        start = MPI_Wtime();
        const auto det_b = determinant(b, MPI_COMM_WORLD);
	if (rank == 0)
        print_determinant("det(B)", det_b, seconds(start));

	if (rank == 0)
        std::cout << "Computing det(C)\n" << std::flush;
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
