#ifndef MODEL_HPP
#define MODEL_HPP
#include "mesh3D.hpp"
#include "tridiag.hpp"
#include <cstddef>
#include <iostream>
#include <vector>

using dvector = std::vector<double>;
/**
 * @brief Abstract base Model class representing the system to solve
 *
 */
class Model {
public:
  /**
   * @brief  Construct a new Model object
   *
   * @param Nx Number of steps in the space discretization
   * @param x0 Start of the discretization
   * @param xend End ot the discretization
   */
  Model(std::size_t Nx, double x0, double xend);
  Model(const Model &m);
  Model &operator=(const Model &m);
  virtual ~Model(){};

  /**
   * @brief Factorize the matrix of the model in a LU format via the Tridiag
   * class
   *
   */
  virtual void setLU() = 0;
  /**
   * @brief Sets the RHS of the model
   *
   */
  virtual void setB(dvector, double, std::size_t, double) = 0;
  /**
   * @brief Sets the initial conditions of the model
   *
   * @param u0 Constant initial condition
   */
  virtual void setU0(double u0);
  virtual void setXend(double d); ///< End of discretization
  /**
   * @brief Set the Mesh object used on this model
   *
   * @param mesh A Mesh object compatible with the model
   */
  virtual void setMesh(Mesh3D &mesh) { M_mesh = &mesh; }
  /**
   * @brief Sets the time step
   *
   * @param d Time step
   */
  void setDt(double d); ///< Time step value

  virtual double b(std::size_t i) const;
  /**
   * @brief Solve the exact solution of the model at a point x
   *
   * @return f(x)
   */
  virtual double solveExact(double);
  std::size_t Nx() const;        ///< Numbers of space steps
  double x0() const;             ///< Start of discretization
  double xend() const;           ///< End of discretization
  double x(std::size_t i) const; ///< Value in vector discretization
  double u(std::size_t i) const; ///< Value in vector solution
  double dx() const;             ///< Space step value
  Tridiag<double> S();           ///< Matrix of the model
  dvector &u();                  ///< Returns whole solution vector
  dvector &b();                  ///< Returns whole RHS vector

protected:
  std::size_t M_Nx;    ///< Number of space steps
  double M_x0;         ///< Start of the discretization
  double M_xend;       ///< End of the discretization
  dvector M_x;         ///< Vector of the discretization
  dvector M_u;         ///< Vector of initial contition
  Tridiag<double> M_S; ///< Matrix of the model
  dvector M_b;         ///< RHS of the model
  double M_dt;         ///< Time step
  Mesh3D *M_mesh;      ///< Pointer on a mesh for 3d visualization
};

#endif
