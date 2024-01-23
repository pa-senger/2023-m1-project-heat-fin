#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

/**
 * @brief A square Tridiagonal matrix
 *
 */
//* With a Tridiag LU matrix Ax=b can be solved with a complexity in O(n)
template <class T> class Tridiag {
public:
  /**
   * @brief Construct a new Tridiag object of size nxn
   *
   * @param n Size of the square matrix
   */
  Tridiag(std::size_t n = 1);
  /**
   * @brief Construct a new Tridiag object by copy
   *
   * @param t Tridiag object to copy
   */
  Tridiag(const Tridiag &t);
  /**
   * @brief Construct a new Tridiag object
   *
   * @param n Size
   * @param l Same value to put in the lower diag
   * @param d Same value to put in the main diag
   * @param u Same value to put in the upper diag
   */
  Tridiag(std::size_t n, T l, T d, T u);
  Tridiag(std::vector<T> l, std::vector<T> d, std::vector<T> u);
  /**
   * @brief Assign operator
   *
   * @param t Tridiag object to copy
   * @return Tridiag& object
   */
  Tridiag &operator=(const Tridiag &t);
  ~Tridiag(){};

  //* Getters & Setters
  std::size_t size() const { return M_n; }
  T upper(int i) const { return M_upper[i]; }
  T diag(int i) const { return M_diag[i]; }
  T lower(int i) const { return M_lower[i]; }
  T &upper(int i) { return M_upper[i]; }
  T &diag(int i) { return M_diag[i]; }
  T &lower(int i) { return M_lower[i]; }
  T operator()(int i, int j) const;
  T &operator()(int i, int j);

  //* Arithmetic Operators & others methods
  /**
   * @brief Do a LU factorization of the Tridiag matrix
   * @details A=LU, since A is tridiagonal L is bidiagonal lower with ones on
   * the diag and U is bidiagonal superior. There is no point storing the ones
   * of L, so we store the whole LU factorization in a single Tridiag object
   * where :
   * lower = the lower part of L
   * diag = main diag of U
   * upper = upper diag of U
   * @return Tridiag object
   */
  Tridiag factorize() const;

  Tridiag &operator+=(const Tridiag &rhs);
  Tridiag &operator-=(const Tridiag &rhs);
  std::vector<T> operator*(const std::vector<T> &rhs) const;
  void clear(); ///< Clear the matrix vectors
  //* The following operators are defined outside the class :
  // Tridiag<T> operator+(const Tridiag<T> &lhs, const Tridiag<T> &rhs)
  // Tridiag<T> operator-(const Tridiag<T> &lhs, const Tridiag<T> &rhs
  // Tridiag<T> operator*(const double lhs, const Tridiag<T> &rhs)
  // Tridiag<T> operator*(const Tridiag<T> &lhs, const double rhs)
  // Tridiag<T> operator/(const Tridiag<T> &lhs, const double rhs)
  //* A method for solving LUx=b is defined outside the class:
  // void solveLU(std::vector<T> &x, const Tridiag<T> &LU, const std::vector<T>
  // &b)

  //* Ostream outsite of the class
  // template <class U>
  // friend std::ostream &operator<<(std::ostream &os, const Tridiag<U> &t);

private:
  std::size_t M_n;        ///< Dimention of the square matrix
  std::vector<T> M_upper; ///< Vector of the upper diagonal part
  std::vector<T> M_diag;  ///< Vector of the main diagonal part
  std::vector<T> M_lower; ///< Vector of the lower diagonal part
};

template <class T>
Tridiag<T>::Tridiag(std::size_t n)
    : M_n(n), M_upper(n - 1, 0.0), M_diag(n, 0.0), M_lower(n - 1, 0.0) {
  // assert(n > 1) : upper and lower are of size n-1 ....
}

template <class T>
Tridiag<T>::Tridiag(const Tridiag &t)
    : M_n(t.M_n), M_upper(t.M_upper), M_diag(t.M_diag), M_lower(t.M_lower) {}

template <class T>
Tridiag<T>::Tridiag(std::size_t n, T l, T d, T u)
    : M_n(n), M_upper(n - 1, u), M_diag(n, d), M_lower(n - 1, l) {}

template <class T>
Tridiag<T>::Tridiag(std::vector<T> l, std::vector<T> d, std::vector<T> u)
    : M_n(d.size()), M_upper(), M_diag(), M_lower() {
  assert(l.size() == M_n - 1);
  assert(u.size() == M_n - 1);

  M_lower = l;
  M_diag = d;
  M_upper = u;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Tridiag<T> &t) {
  for (std::size_t i = 0; i < t.size(); ++i) {
    for (std::size_t j = 0; j < t.size(); ++j) {
      if (i == j)
        os << t.diag(i) << "\t";
      else if (i == j + 1)
        os << t.lower(i - 1) << "\t";
      else if (i == j - 1)
        os << t.upper(i) << "\t";
      else
        os << ". \t";
    }
    os << "\n";
  }
  return os;
}

template <class T> Tridiag<T> &Tridiag<T>::operator=(const Tridiag &t) {
  if (this != &t) {
    M_n = t.M_n;
    M_upper = t.M_upper;
    M_diag = t.M_diag;
    M_lower = t.M_lower;
  }
  return *this;
}

template <class T> Tridiag<T> Tridiag<T>::factorize() const {
  Tridiag<T> LU(M_n);

  LU.diag(0) = M_diag[0];
  for (std::size_t i = 1; i < M_n; ++i) {
    LU.lower(i - 1) = M_lower[i - 1] / LU.diag(i - 1);
    LU.upper(i - 1) = M_upper[i - 1];
    LU.diag(i) = M_diag[i] - LU.lower(i - 1) * LU.upper(i - 1);
  }
  return LU;
}

template <class T> Tridiag<T> &Tridiag<T>::operator-=(const Tridiag &rhs) {
  assert(M_n == rhs.size());

  for (std::size_t i = 0; i < M_n - 1; ++i) {
    M_lower[i] -= rhs.lower(i);
    M_upper[i] -= rhs.upper(i);
    M_diag[i] -= rhs.diag(i);
  }
  M_diag[M_n - 1] -= rhs.diag(M_n - 1);
  return *this;
}

template <class T> Tridiag<T> &Tridiag<T>::operator+=(const Tridiag &rhs) {
  assert(M_n == rhs.size());

  for (std::size_t i = 0; i < M_n - 1; ++i) {
    M_lower[i] += rhs.lower(i);
    M_upper[i] += rhs.upper(i);
    M_diag[i] += rhs.diag(i);
  }
  M_diag[M_n - 1] += rhs.diag(M_n - 1);
  return *this;
}

template <class T>
inline Tridiag<T> operator+(const Tridiag<T> &lhs, const Tridiag<T> &rhs) {
  assert(lhs.size() == rhs.size());

  Tridiag<T> res = lhs;
  res += rhs;
  return res;
}

template <class T>
inline Tridiag<T> operator-(const Tridiag<T> &lhs, const Tridiag<T> &rhs) {
  assert(lhs.size() == rhs.size());

  Tridiag<T> res = lhs;
  res -= rhs;
  return res;
}

template <class T>
inline Tridiag<T> operator*(const double lhs, const Tridiag<T> &rhs) {
  Tridiag<T> res = rhs;
  for (std::size_t i = 0; i < rhs.size() - 1; ++i) {
    res.lower(i) = rhs.lower(i) * lhs;
    res.upper(i) = rhs.upper(i) * lhs;
    res.diag(i) = rhs.diag(i) * lhs;
  }
  res.diag(rhs.size() - 1) = rhs.diag(rhs.size() - 1) * lhs;

  return res;
}

template <class T>
inline Tridiag<T> operator*(const Tridiag<T> &lhs, const double rhs) {
  Tridiag<T> res;
  res = rhs * lhs;
  return res;
}

template <class T>
inline Tridiag<T> operator/(const Tridiag<T> &lhs, const double rhs) {
  Tridiag<T> res = lhs;
  for (std::size_t i = 0; i < lhs.size() - 1; ++i) {
    res.lower(i) = lhs.lower(i) / rhs;
    res.upper(i) = lhs.upper(i) / rhs;
    res.diag(i) = lhs.diag(i) / rhs;
  }
  res.diag(lhs.size() - 1) = lhs.diag(lhs.size() - 1) / rhs;

  return res;
}

template <class T>
void solveL(std::vector<T> &y, const Tridiag<T> &LU, const std::vector<T> &b) {
  assert(y.size() == LU.size());
  assert(b.size() == LU.size());

  y[0] = b[0];
  for (std::size_t i = 1; i < LU.size(); ++i)
    y[i] = b[i] - LU.lower(i - 1) * y[i - 1];
}

template <class T>
void solveU(std::vector<T> &x, const Tridiag<T> &LU, const std::vector<T> &y) {
  assert(x.size() == LU.size());
  assert(y.size() == LU.size());

  x[LU.size() - 1] = y[LU.size() - 1] / LU.diag(LU.size() - 1);
  for (int i = LU.size() - 2; i >= 0; --i) {
    if (std::fabs(LU.diag(i)) < 1e-8)
      std::cerr << "Error: pivot is too close to zero" << std::endl;
    x[i] = (y[i] - LU.upper(i) * x[i + 1]) / LU.diag(i);
  }
}
/**
 * @brief Solve a LUx=b system in O(n) time complexity
 *
 * @tparam T Typename
 * @param x Vector of solution to fill
 * @param LU Tridiag matrix factorized in LU format
 * @param b RHS of the system
 */
template <class T>
void solveLU(std::vector<T> &x, const Tridiag<T> &LU, const std::vector<T> &b) {
  assert(x.size() == LU.size());
  assert(b.size() == LU.size());

  std::vector<T> y(LU.size());
  solveL(y, LU, b);
  solveU(x, LU, y);
}

template <class T> void Tridiag<T>::clear() {
  M_upper.clear();
  M_diag.clear();
  M_lower.clear();
  M_upper.resize(M_n - 1);
  M_diag.resize(M_n);
  M_lower.resize(M_n - 1);
}

template <class T> T Tridiag<T>::operator()(int i, int j) const {
  if (i == j)
    return diag(i);
  if (i == j + 1)
    return lower(i - 1);
  if (i == j - 1)
    return upper(i);

  return T(); // if not in the band
}
template <class T> T &Tridiag<T>::operator()(int i, int j) {
  if (i == j)
    return M_diag[i];
  if (i == j + 1)
    return M_lower[i - 1];
  if (i == j - 1)
    return M_upper[i];
  std::cerr << "Error: the coef (" + std::to_string(i) + ", " +
                   std::to_string(j) + ") is not in the band of the matrix "
            << std::endl;
  exit(1);
  return M_diag[0]; // to avoid warning at compilation
}

template <class T>
std::vector<T> Tridiag<T>::operator*(const std::vector<T> &rhs) const {
  std::vector<T> res(M_n, 0);

  // Diagonal multiplication
  for (std::size_t i = 0; i < M_n; ++i)
    res[i] += M_diag[i] * rhs[i];

  // Upper diagonal multiplication
  for (std::size_t i = 0; i < M_n - 1; ++i)
    res[i] += M_upper[i] * rhs[i + 1];

  // Lower diagonal multiplication
  for (std::size_t i = 1; i < M_n; ++i)
    res[i] += M_lower[i - 1] * rhs[i - 1];

  return res;
};

#endif
