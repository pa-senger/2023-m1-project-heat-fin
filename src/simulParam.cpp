#include "../include/simulParam.hpp"
#include <fstream>
#include <iostream>

// Function 'operator<<' defined in a header file; function definitions in
// header files can lead to ODR violations

std::ostream &operator<<(std::ostream &os, const SimulParam &params) {
  os << "Lx " << params.Lx << " Ly " << params.Ly << " Lz " << params.Lz
     << "\n";
  os << "Nx " << params.Nx << "\n";

  os << "Phi " << params.Phi << " hc " << params.hc << " Te " << params.Te
     << "\n";
  os << "rho " << params.rho << " kappa " << params.kappa << "\n";

  os << "stationary " << params.stationary << " cycling " << params.cycling
     << "\n";
  os << "fan " << params.fan << " cooling " << params.cooling << "\n";

  os << "tFinal " << params.tFinal << " Nt " << params.Nt << "\n";

  os << "Mx " << params.Mx << " My " << params.My << " Mz " << params.Mz
     << "\n";
  os << "doPlots " << params.doPlots << " do3D " << params.do3D << "\n";
  os << "solName " << params.solName << "\n";
  return os;
}
std::istream &operator>>(std::istream &is, SimulParam &params) {
  ///< Assumes the input format is always consistent
  std::string str;
  double val;
  is >> str >> val;
  params.Lx = val * 1e-3;

  is >> str >> val;
  params.Ly = val * 1e-3;

  is >> str >> val;
  params.Lz = val * 1e-3;

  is >> str >> params.Nx;

  is >> str >> val;
  params.Phi = val * 1e6;

  is >> str >> val;
  params.hc = val * 1e6;

  is >> str >> params.Te;
  is >> str >> params.rho;
  is >> str >> params.kappa;
  is >> str >> params.stationary;
  is >> str >> params.cycling;
  is >> str >> params.fan;
  is >> str >> params.cooling;
  is >> str >> params.tFinal;
  is >> str >> params.Nt;
  is >> str >> params.Mx >> str >> params.My >> str >> params.Mz;
  is >> str >> params.doPlots >> str >> params.do3D;
  is >> str >> params.solName;
  return is;
}
