add_executable(algebra_determinants main.cpp)
target_compile_features(algebra_determinants PRIVATE cxx_std_20)
target_compile_definitions(algebra_determinants PRIVATE
    ALGEBRA_ASSIGNMENT_SIZE=${ALGEBRA_ASSIGNMENT_SIZE}
)
target_link_libraries(algebra_determinants PRIVATE
    algebra::finite_fields
    algebra::matrices
)
find_package(MPI COMPONENTS CXX)
if(MPI_CXX_FOUND)
    add_executable(mpi_efficient mpi_efficient.cpp)
    target_compile_features(mpi_efficient PRIVATE cxx_std_20)
    target_compile_definitions(mpi_efficient PRIVATE ALGEBRA_ASSIGNMENT_SIZE=${ALGEBRA_ASSIGNMENT_SIZE})
    target_link_libraries(mpi_efficient PRIVATE algebra::finite_fields algebra::matrices MPI::MPI_CXX)

    add_executable(mpi_inefficient mpi_inefficient.cpp)
    target_compile_features(mpi_inefficient PRIVATE cxx_std_20)
    target_compile_definitions(mpi_inefficient PRIVATE ALGEBRA_ASSIGNMENT_SIZE=${ALGEBRA_ASSIGNMENT_SIZE})
    target_link_libraries(mpi_inefficient PRIVATE algebra::finite_fields algebra::matrices MPI::MPI_CXX)
else()
    message(STATUS "MPI not found: skipping the MPI determinant programs")
endif()