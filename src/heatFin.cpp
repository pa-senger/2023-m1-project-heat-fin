#include "../include/heatFin.hpp"
#include <cstddef>
#include <ostream>

HeatFin::HeatFin(double nx, double rho, double c, double k, double Te,
                 double phi, double hc, double Lx, double Ly, double Lz)
    : Model(nx, 0, Lx), M_rho(rho), M_C(c), M_k(k), M_Te(Te), M_phi(phi),
      M_hc(hc), M_Ly(Ly), M_Lz(Lz), M_cooler(0), M_flux(true), M_cooling(false),
      M_cyclingStatus(false), M_stationnary(false), M_staticInterpolation() {}

HeatFin::HeatFin(const HeatFin &f)
    : Model(f), M_rho(f.M_rho), M_C(f.M_C), M_k(f.M_k), M_Te(f.M_Te),
      M_phi(f.M_phi), M_hc(f.M_hc), M_Ly(f.M_Ly), M_Lz(f.M_Lz),
      M_cooler(f.M_cooler), M_flux(f.M_flux), M_cooling(f.M_cooling),
      M_cyclingStatus(f.M_cyclingStatus), M_stationnary(f.M_stationnary),
      M_staticInterpolation(f.M_staticInterpolation) {}

HeatFin &HeatFin::operator=(const HeatFin &f) {
  if (this != &f) {
    Model::operator=(f);
    M_rho = f.M_rho;
    M_C = f.M_C;
    M_k = f.M_k;
    M_Te = f.M_Te;
    M_phi = f.M_phi;
    M_hc = f.M_hc;
    M_Ly = f.M_Ly;
    M_Lz = f.M_Lz;
    M_flux = f.M_flux;
    M_cooler = f.M_cooler;
    M_cooling = f.M_cooling;
    M_cyclingStatus = f.M_cyclingStatus;
    M_stationnary = f.M_stationnary;
  }
  return *this;
}

void HeatFin::setLU() {
  Tridiag<double> P(M_Nx + 1, -1, 2, -1);
  M_S = P;
  // M_S.upper(0) = -dx();
  // M_S.lower(M_Nx - 1) = -dx();
  // M_S.diag(0) = dx();
  // M_S.diag(M_Nx) = dx();
  M_S(0, 0) = dx();
  M_S(M_Nx, M_Nx) = dx();
  M_S(0, 1) = -dx();
  M_S(M_Nx, M_Nx - 1) = -dx();

  Tridiag<double> I(M_Nx + 1, 0, 1, 0);
  I(0, 0) = 0;
  I(M_Nx, M_Nx) = 0;

  Tridiag<double> LU(M_Nx + 1);
  if (M_stationnary) {
    LU = (M_k / pow(dx(), 2)) * M_S + (M_hc * p() / s()) * I;
  } else {
    LU = (M_rho * M_C / M_dt + M_hc * p() / s()) * I +
         (M_k / pow(dx(), 2) * M_S);
  }
  M_S = LU.factorize();
}

void HeatFin::setB(dvector u, double t, std::size_t nt, double tf) {
  M_b.clear();
  M_b.resize(M_Nx + 1);
  if (M_stationnary) {
    for (std::size_t i = 1; i < M_Nx; ++i)
      M_b[i] = M_Te * M_hc * p() / s();

  } else {
    double cycle = 30;
    for (std::size_t i = 1; i < M_Nx; ++i)
      M_b[i] = M_rho * M_C * u[i] / M_dt + M_Te * M_hc * p() / s();

    if (M_cyclingStatus) {
      int currentCycle = t / (cycle * nt / tf);
      (currentCycle % 2 == 0) ? fluxON() : fluxOFF();
    }
  }
  M_b[M_Nx] = 0;
  M_flux ? M_b[0] = M_phi : M_b[0] = 0;
  cooling() ? CoolingON() : CoolingOFF();
}

double HeatFin::solveExact(double x) {
  double a = (M_hc * p()) / (M_k * s());
  double num = M_phi * cosh(sqrt(a) * M_xend) * cosh(sqrt(a) * (M_xend - x));
  double denom =
      M_k * sqrt(a) * sinh(sqrt(a) * M_xend) * cosh(sqrt(a) * M_xend);
  return M_Te + num / denom;
}

void HeatFin::fanON() { M_hc = 200; }

void HeatFin::fanOFF() { M_hc = 10; }

bool HeatFin::flux() const { return M_flux; }

bool HeatFin::cooling() const { return M_cooling; }

void HeatFin::fluxON() {
  M_flux = true;
  M_b[0] = M_phi;
}

void HeatFin::fluxOFF() {
  M_flux = false;
  M_b[0] = 0;
}

void HeatFin::CoolingON() {
  M_cooling = true;
  M_phi += M_cooler;
}

void HeatFin::CoolingOFF() { M_cooling = false; }

double HeatFin::p() const { return 2 * (M_Ly + M_Lz); }

double HeatFin::s() const { return M_Lz * M_Ly; }

void HeatFin::setHc(double power) { M_hc = power; }

double HeatFin::hc() const { return M_hc; }

void HeatFin::setCooler(double temp) { M_cooler = temp; }

void HeatFin::setPhi(double d) { M_phi = d; }

void HeatFin::setLx(double d) { setXend(d); }

void HeatFin::setCycling(bool b) { M_cyclingStatus = b; }

void HeatFin::setStationnary(bool b) { M_stationnary = b; }

void HeatFin::saveStaticInterpolation(const std::string filename) {
  if (M_mesh == nullptr) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }
  std::filesystem::path filePath(filename);

  // Ensure the directory exists
  std::filesystem::create_directories(filePath.parent_path());
  std::ofstream ofile(filePath, std::ios::out | std::ios::trunc);

  if (ofile) {
    ofile << "# vtk DataFile Version 2.0\n";
    ofile << "vtk output\n";
    ofile << "ASCII\n";
    ofile << "DATASET STRUCTURED_GRID\n";
    ofile << "DIMENSIONS " << M_mesh->Mx() << " " << M_mesh->My() << " "
          << M_mesh->Mz() << "\n";
    ofile << "POINTS " << M_mesh->Mx() * M_mesh->My() * M_mesh->Mz()
          << " float\n";
    for (std::size_t k = 0; k < M_mesh->Mz(); ++k) {
      for (std::size_t j = 0; j < M_mesh->My(); ++j) {
        for (std::size_t i = 0; i < M_mesh->Mx(); ++i) {
          ofile << i << " " << j << " " << k << "\n";
        }
      }
    }
    ofile << "POINT_DATA " << M_mesh->Mx() * M_mesh->My() * M_mesh->Mz()
          << " \n";
    ofile << "FIELD FieldData 1\n";
    ofile << "sol1 1 " << M_mesh->Mx() * M_mesh->My() * M_mesh->Mz()
          << " float\n";
    for (std::size_t k = 0; k < M_mesh->Mz(); ++k) {
      for (std::size_t j = 0; j < M_mesh->My(); ++j) {
        for (std::size_t i = 0; i < M_mesh->Mx(); ++i) {
          ofile << M_staticInterpolation[i] << "\n";
        }
      }
    }
    ofile << std::flush;
  } else {
    std::cerr << "Error: cant open " + filename + " !" << std::endl;
    exit(1);
  }
  // std::string instruction = "paraview " + filename;
}

void HeatFin::interpolate(dvector &U) {
  if (M_mesh == nullptr) {
    std::cerr << "Error: You neet to set a Mesh first" << std::endl;
    exit(1);
  }
  M_staticInterpolation.clear();
  for (std::size_t i = 0; i <= M_mesh->Mx(); ++i) {
    std::size_t k = 0;
    while (M_x[k] < M_mesh->x(i)) {
      ++k;
    }
    double b = U[0];
    double a = (U[k] + U[k + 1] - 2 * b) / (M_x[k] + M_x[k + 1]);
    M_staticInterpolation.push_back(a * M_mesh->x(i) + b);
  }
}

void HeatFin::interpolate(std::vector<dvector> &T, std::size_t Nt,
                          const std::string solName) {
  for (std::size_t t = 0; t <= Nt; t++) {
    interpolate(T[t]);
    saveStaticInterpolation("../data/3d/" + solName + "." + std::to_string(t) +
                            ".vtk");
  }
}