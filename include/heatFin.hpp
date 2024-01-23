#ifndef HEAT_FIN_STATIC_HPP
#define HEAT_FIN_STATIC_HPP
#include "model.hpp"
#include "solverTime.hpp"
#include "tridiag.hpp"
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

/**
 * @brief A Model of a single heatfin of a CPU radiator
 *
 */
class HeatFin : public Model {
public:
  /**
   * @brief Construct a new Heat Fin object
   *
   * @param nx Number of space steps in the discretization
   * @param rho Density of the material
   * @param c Heat at constant pressure
   * @param k Thermal condictivity
   * @param Te External tempeature
   * @param phi Heat flux power
   * @param hc Heat surfacic transfer coefficient
   * @param Lx Size of the fin on the x-axis
   * @param Ly Size of the fin on the y-axis
   * @param Lz Size of the fin on the z-axis
   */
  HeatFin(double nx, double rho, double c, double k, double Te, double phi,
          double hc, double Lx, double Ly, double Lz);
  HeatFin(const HeatFin &f);
  HeatFin &operator=(const HeatFin &f);
  ~HeatFin() {}

  void setLU() override;
  void setB(dvector u, double t, std::size_t nt, double tf) override;
  void fanON();
  void fanOFF();
  void fluxON();
  void fluxOFF();
  void CoolingON();
  void CoolingOFF();
  void setHc(double power);    ///< Fan power
  void setCooler(double temp); ///< Cooling power
  void setPhi(double d);       ///< Heat flux power
  void setLx(double d);        ///< Same as setXend() from base class
  void setCycling(bool b);     ///< Cycling status ON or OFF
  void setStationnary(bool b); ///< Dynamic or stationnary model

  /**
   * @brief Interpolate the stationnary 1d model on a mesh for 3d visualisation
   *
   * @param U Vector of 1d solution to sets on the mesh
   */
  void interpolate(dvector &U);

  /**
   * @brief Interpolate the dynamic solutions on the mesh
   *
   * @param T Vector of vectors of solutions in time
   * @param Nt Number of time steps
   * @param solName Name to save in /data/3d
   */
  void interpolate(std::vector<dvector> &T, std::size_t Nt,
                   const std::string solName);

  /**
   * @brief Saves the 3d visualisation file of the stationnary model
   *
   * @param filename Name of the file, including path and extension
   */
  void saveStaticInterpolation(const std::string filename);

  double solveExact(double x) override;
  bool flux() const;    ///< Current status of the flux, cycling or not
  bool cooling() const; ///< Current status of the cooling
  double p() const;     ///< Perimeter
  double s() const;     ///< Surface
  double hc() const;    ///< Fan power

protected:
  double M_rho; ///< Density
  double M_C;   ///< Heat at constant pressure
  double M_k;   ///< Thermal conductivity
  double M_Te;  ///< External Temperature
  double M_phi; ///< Heat flux power
  double M_hc;  ///< Heat surfacic transfer coefficient
  double M_Ly;  ///< Size of the fin on the y-axis
  double M_Lz;  ///< Size of the fin on the z-axis
  ///< M_Lx = Model::M_xend
  double M_cooler;      ///< Cooling power effect (so <0)
  bool M_flux;          ///< Heat flux status
  bool M_cooling;       ///< Cooling status, why did I do this ...?
  bool M_cyclingStatus; ///< Status of the heat flux
  bool M_stationnary;   ///< Stationnary or Dynamic model
  dvector M_staticInterpolation;
  // Cooling need to be changed to accuratly describe what happens IRL
};

#endif
