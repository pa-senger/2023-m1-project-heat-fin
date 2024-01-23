#include "../include/model.hpp"
#include "../include/tridiag.hpp"
#include <cstddef>

Model::Model(std::size_t Nx, double x0, double xend)
    : M_Nx(Nx), M_x0(x0), M_xend(xend), M_x(Nx + 1), M_u(Nx + 1),
      M_S(Nx + 1, -1, 2, -1), M_b(Nx + 1), M_dt(0), M_mesh(nullptr) {
  for (std::size_t i = 0; i <= M_Nx; ++i)
    M_x[i] = M_x0 + i * dx();
}

Model::Model(const Model &m)
    : M_Nx(m.M_Nx), M_x0(m.M_x0), M_xend(m.M_xend), M_x(m.M_x), M_u(m.M_u),
      M_S(m.M_S), M_b(m.M_b), M_dt(m.M_dt), M_mesh(m.M_mesh) {}

Model &Model::operator=(const Model &m) {
  if (this != &m) {
    M_Nx = m.M_Nx;
    M_x0 = m.M_x0;
    M_xend = m.M_xend;
    M_x = m.M_x;
    M_u = m.M_u;
    M_S = m.M_S;
    M_b = m.M_b;
    M_dt = m.M_dt;
    M_mesh = m.M_mesh;
  }
  return *this;
}

double Model::b(std::size_t i) const { return M_b[i]; }

void Model::setU0(double d) {
  for (std::size_t i = 0; i < M_Nx + 1; ++i)
    M_u[i] = d;
}

double Model::solveExact(double) { return 0; }

void Model::setXend(double d) {
  M_xend = d;
  for (std::size_t i = 0; i < M_Nx + 1; ++i)
    M_x[i] = M_x0 + i * dx();
}

void Model::setDt(double d) { M_dt = d; }

std::size_t Model::Nx() const { return M_Nx; }

double Model::x0() const { return M_x0; }

double Model::xend() const { return M_xend; }

double Model::x(std::size_t i) const { return M_x[i]; }

double Model::u(std::size_t i) const { return M_u[i]; }

double Model::dx() const { return (M_xend - M_x0) / M_Nx; }

Tridiag<double> Model::S() { return M_S; }

dvector &Model::u() { return M_u; }

dvector &Model::b() { return M_b; }
