#include "../include/tridiag.hpp"
#include <cmath>
#include <vector>

#define TOL 1e-4

int main() {
  // Create a Tridiag object P like Poisson with the specified values
  Tridiag<double> P(3, -1, 2, -1);

  // Factorize the Tridiag object P
  Tridiag<double> LU = P.factorize();
  std::cout << "LU=\n" << LU << std::endl;

  // Assertions to check the factorization
  assert(std::fabs(LU.diag(0) - 2) < TOL);
  assert(std::fabs(LU.diag(1) - 1.5) < TOL);
  assert(std::fabs(LU.diag(2) - 1.33333) < TOL);
  assert(std::fabs(LU.upper(0) + 1) < TOL);
  assert(std::fabs(LU.upper(1) + 1) < TOL);
  assert(std::fabs(LU.lower(0) + 0.5) < TOL);
  assert(std::fabs(LU.lower(1) + 0.666667) < TOL);

  // Create input vectors b and x
  std::vector<double> b = {1., 1., 2.};
  std::vector<double> x(3);

  // Solve for x using the LU factorization
  solveLU(x, LU, b);

  // Assertions to check the solution
  assert(std::fabs(x[0] - 1.7499999999999998) < TOL);
  assert(std::fabs(x[1] - 2.4999999999999996) < TOL);
  assert(std::fabs(x[2] - 2.2499999999999996) < TOL);

  // Tests Tridiag * vector product
  std::vector<double> l = {6, 1.5, 4};
  std::vector<double> d = {1, 4, 2, 8};
  std::vector<double> u = {2, 5, 9};

  Tridiag<double> A(l, d, u);

  // Initialize vector b
  std::vector<double> b2 = {7.9, 2.5, 3.3, 11};

  // Expected result
  std::vector<double> expected_result = {12.9, 73.9, 109.35, 101.2};

  // Perform matrix-vector multiplication
  std::vector<double> result = A * b2;

  // Test using assert
  for (std::size_t i = 0; i < result.size(); ++i) {
    assert(std::fabs(result[i] - expected_result[i]) < TOL);
  }

  return 0;
}
