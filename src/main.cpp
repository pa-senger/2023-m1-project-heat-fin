#include "../include/heatFin.hpp"
#include "../include/mesh3D.hpp"
#include "../include/model.hpp"
#include "../include/simulParam.hpp"
#include "../include/solverTime.hpp"
#include "../include/tridiag.hpp"
#include <chrono>
#include <cmath>
#include <cstddef>
#include <vector>

//* To compile : cd in a fresh build directory then do :
// cmake --preset release ..
// make
//* To run the program do  : ./run <config_file> ; for example :
// ./run ../simul.cfg
//* To view some graphs and testings (you need pandas and matplotlib) do :
// make test

using dvector = std::vector<double>;

int main(int argc, char **argv) {
  // Start a chrono
  std ::chrono ::time_point<std ::chrono ::system_clock> start, end;
  start = std ::chrono ::system_clock ::now();

  //* Set configurations for the simulation

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
    return 1;
  }

  // Read simulation parameters from the cfg file
  std::ifstream configFile(argv[1]);
  if (!configFile.is_open()) {
    std::cerr << "Error: Unable to open the configuration file " << argv[1]
              << std::endl;
    return 1;
  }
  // Put parameters into the Struct
  SimulParam params;
  configFile >> params;
  configFile.close();
  bool plot = params.doPlots;

  // Print parameters of the simultation into the console
  std::cout << params << std::endl;

  //* Physical paramameters
  double Cp = 940; // Heat at constant pressure

  //* Run simulation
  HeatFin f(params.Nx, params.rho, Cp, params.kappa, params.Te, params.Phi,
            params.hc, params.Lx, params.Ly, params.Lz);
  f.setStationnary(params.stationary);
  f.setCycling(params.cycling);
  params.fan ? f.fanON() : f.fanOFF();
  params.cooling ? f.CoolingON() : f.CoolingOFF();

  SolverTime s1(params.Nt);
  s1.setModel(f);
  Mesh3D mesh(params.Lx, params.Ly, params.Lx, params.Mx, params.My, params.Mz);
  f.setMesh(mesh);

  //* Stationnary Model
  if (params.stationary) {
    s1.solveStatic();
    s1.saveStatic("../data/2d/" + params.solName + ".csv", plot);

    if (params.do3D) {
      f.interpolate(s1.T().back());
      f.saveStaticInterpolation("../data/3d/" + params.solName + ".vtk");
    }
  } else {
    //* Dynamic Model
    dvector times;
    dvector points(3);
    double step = 15;
    for (double t = step; t <= params.tFinal; t += 1.5 * step)
      times.push_back(t);
    points = {0, params.Lx / 2, params.Lx};

    f.setU0(params.Te); // Initial Conditions
    s1.solve(params.tFinal);
    // ~=Constant solutions at certain times are expected in a cycling model
    s1.saveAtTimes("../data/2d/times_" + params.solName + ".csv", times, plot);
    s1.saveAtPoints("../data/2d/points_" + params.solName + ".csv", points,
                    plot);

    if (params.do3D)
      f.interpolate(s1.T(), params.Nt, params.solName);
  }

  // End chrono
  end = std::chrono::system_clock::now();

  // Print simulation time
  double elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  double elapsed_s =
      std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  std::cout << "The simulation took: " << elapsed_ms << "ms (approx "
            << elapsed_s << "s)\n";

  return 0;
}
