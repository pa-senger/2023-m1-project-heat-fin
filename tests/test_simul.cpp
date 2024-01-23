#include "../include/heatFin.hpp"
#include "../include/model.hpp"
#include "../include/simulParam.hpp"
#include "../include/solverTime.hpp"
#include "../include/tridiag.hpp"
#include <cmath>
#include <cstddef>
#include <vector>

int main() {

  //* Physical paramameters
  double rho = 2700;
  double Cp = 940;
  double k = 164;
  double Te = 20;
  double Phi = 1.25 * 1e5;
  double hc = 200;
  double Lx = 0.04;
  double Ly = 0.004;
  double Lz = 0.05;

  std::size_t Nx = 10'000;
  std::size_t Nt = 600;
  double tFinal = 300;

  //* Stationnary Model
  HeatFin f(Nx, rho, Cp, k, Te, Phi, hc, Lx, Ly, Lz);
  f.setStationnary(1);
  SolverTime s1(Nt);
  s1.setModel(f);
  s1.solveStatic();
  s1.saveStatic("../data/test/static.csv", true);
  f.setLx(2 * Lx);
  s1.solveStatic();
  s1.saveStatic("../data/test/static.csv", true);
  f.setLx(Lx);

  //* Dynamic model

  f.setStationnary(0);

  f.setU0(Te); // Initial Conditions
  s1.solve(tFinal);
  std::vector<double> times(6);
  std::vector<double> points(3);
  times = {15, 30, 65, 90, 150, 210};
  points = {0, Lx / 2, Lx};
  s1.saveAtTimes("../data/test/times.csv", times, true);
  s1.saveAtTimes("../data/test/times.csv", true, 0, 10, 15, 150, 200);
  f.setCycling(1);
  s1.solve();
  s1.saveAtPoints("../data/test/points.csv", true, 0.02, 0.03);
  s1.saveAtPoints("../data/test/points.csv", points, true);
  SolverTime s2(s1);
  s2.solve();
  s1.saveAtPoints("../data/test/points.csv", points, true);
  HeatFin f2 = f;
  s2.setModel(f2);
  f2.fanOFF();
  s2.solve();
  s2.saveAtPoints("../data/test/fanOFF.csv", points, true);
  f2.fanON();
  f2.setCooler(-200);
  f2.CoolingON();
  s2.solve();
  s2.saveAtPoints("../data/test/coolingON.csv", points, true);
  f2.CoolingOFF();
  f2.setPhi(Phi);

  return 0;
}
