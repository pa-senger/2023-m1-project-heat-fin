#ifndef SIMULPARAM_HPP
#define SIMULPARAM_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/**
 * @brief The parameters of the simulation given by a cfg file
 *
 */
struct SimulParam {
  // This contructor exist solely to silence -Weffc++ ...
  SimulParam()
      : Lx(), Ly(), Lz(), Phi(), hc(), Te(), tFinal(), rho(), kappa(), Nx(),
        Nt(), Mx(), My(), Mz(), stationary(), cycling(), doPlots(), fan(),
        cooling(), do3D(), solName() {}

  double Lx, Ly, Lz, Phi, hc, Te, tFinal, rho, kappa;
  std::size_t Nx, Nt, Mx, My, Mz;
  bool stationary, cycling, doPlots, fan, cooling, do3D;
  std::string solName;

  friend std::ostream &operator<<(std::ostream &os, const SimulParam &params);
  friend std::istream &operator>>(std::istream &is, SimulParam &params);
};

#endif