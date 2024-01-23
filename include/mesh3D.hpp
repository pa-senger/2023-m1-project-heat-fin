#ifndef MESH3D_HPP
#define MESH3D_HPP

#include "tridiag.hpp"
#include <cstddef>
#include <vector>

using sizeVector = std::vector<std::size_t>;
using dvector = std::vector<double>;

/**
 * @brief Mesh for a 3d representation of a Model
 * @details This class is not very important yet and mostly useless
 *
 */
class Mesh3D {
public:
  /**
   * @brief Construct a new Mesh3D object
   *
   * @param Lx Size of the model to mesh on the x-axis
   * @param Ly Size of the model to mesh on the y-axis
   * @param Lz Size of the model to mesh on the z-axis
   * @param x Number of space steps on the x-axis of the mesh
   * @param y Number of space steps on the y-axis of the mesh
   * @param z Number of space steps on the z-axis of the mesh
   */
  Mesh3D(double Lx, double Ly, double Lz, int x, int y, int z)
      : M_Lxyz(3), M_Mxyz(3), M_x(x + 1), M_y(y + 1), M_z(z + 1) {
    setLxyz(Lx, Ly, Lz);
    setMxyz(x, y, z);

    for (std::size_t i = 0; i <= Mx(); ++i)
      M_x[i] = i * dx();
    for (std::size_t i = 0; i <= My(); ++i)
      M_y[i] = i * dy();
    for (std::size_t i = 0; i <= Mz(); ++i)
      M_z[i] = i * dz();
  }
  double Lxyz(int i) const { return M_Lxyz[i]; }
  double &Lxyz(int i) { return M_Lxyz[i]; }
  std::size_t Mxyz(int i) const { return M_Mxyz[i]; }
  std::size_t &Mxyz(int i) { return M_Mxyz[i]; }
  double Lx() const { return M_Lxyz[0]; } ///< Size of the model on x-axis
  double Ly() const { return M_Lxyz[1]; } ///< Size of the model on y-axis
  double Lz() const { return M_Lxyz[2]; } ///< Size of the model on z-axis
  double x(int i) const { return M_x[i]; }
  double y(int i) const { return M_y[i]; }
  double z(int i) const { return M_z[i]; }
  double &x(int i) { return M_x[i]; }
  double &y(int i) { return M_y[i]; }
  double &z(int i) { return M_z[i]; }
  std::size_t Mx() const { return M_Mxyz[0]; } ///< Number of space steps
  std::size_t My() const { return M_Mxyz[1]; } ///< Number of space steps
  std::size_t Mz() const { return M_Mxyz[2]; } ///< Number of space steps
  std::size_t &Mx() { return M_Mxyz[0]; }      ///< Number of space steps
  std::size_t &My() { return M_Mxyz[1]; }      ///< Number of space steps
  std::size_t &Mz() { return M_Mxyz[2]; }      ///< Number of space steps
  void setLxyz(double x, double y, double z) {
    M_Lxyz[0] = x;
    M_Lxyz[1] = y;
    M_Lxyz[z] = z;
  }
  void setMxyz(int x, int y, int z) {
    M_Mxyz[0] = x;
    M_Mxyz[1] = y;
    M_Mxyz[2] = z;
  }
  double dx() const { return Lx() / Mx(); }
  double dy() const { return Ly() / My(); }
  double dz() const { return Lz() / Mz(); }
  dvector x_ijk(int i, int j, int k) const {
    dvector res(3);
    res[0] = x(i);
    res[1] = y(j);
    res[2] = z(k);
    return res;
  }

private:
  dvector M_Lxyz;    ///< Dimention of the model to mesh
  sizeVector M_Mxyz; ///< Number of steps in each dimention
  dvector M_x;       ///< Discretization in the x-axis
  dvector M_y;       ///< Discretization in the y-axis
  dvector M_z;       ///< Discretization in the z-axis
};

#endif