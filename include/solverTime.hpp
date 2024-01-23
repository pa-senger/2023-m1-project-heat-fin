#ifndef SOLVERTIME_HPP
#define SOLVERTIME_HPP
#include "model.hpp"
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

using dvector = std::vector<double>;

/**
 * @brief A solver class for the Model class
 *
 */
class SolverTime {
public:
  /**
   * @brief Construct a new Solver Time object
   *
   * @param Nt Number of time steps
   * @param tf Final time
   */
  SolverTime(std::size_t Nt = 0, double tf = 1);
  SolverTime(const SolverTime &s);
  SolverTime &operator=(const SolverTime &s);
  ~SolverTime();

  double dt() const; ///< Get the time step value

  void setTfinal(double t);   ///< Sets the final time value
  void setNt(std::size_t nt); ///< Sets the number of time steps
  void setModel(Model &m);    ///< Sets a pointer to a Model
  void solveStatic();         ///< Solves a stationnary Model
  void solve();               ///< Solves a dynamic Model in time

  /**
   * @brief Solves a dynamic Model in time
   *
   * @param tFinal Final time
   */
  void solve(double tFinal);

  std::vector<dvector> &T() { return M_T; } ///< Gets the vector of solutions

  /**
   * @brief Saves the stationnary 2d solution and do a plot
   * @details The exact solution will also be added in the file
   *
   * @param filename Filename including path and extention
   * @param do_plot Instruction to do the plot or only save data
   */
  void saveStatic(const std::string filename, const bool do_plot = false);

  /**
   * @brief Saves the 2d solution at multiple times
   *
   * @param filename Filename including path and extention
   * @param times Vector of timestamps
   * @param do_plot Instruction to do the plot or only save data
   */
  void saveAtTimes(const std::string filename, dvector times,
                   const bool do_plot = false);

  template <typename... Args>
  void saveAtTimes(const std::string &filename, bool do_plot, Args... args);

  /**
   * @brief Saves the 2d solution at multiple points
   *
   * @param filename Filename including path and extention
   * @param times Vector of points
   * @param do_plot Instruction to do the plot or only save data
   */
  void saveAtPoints(const std::string filename,
                    const std::vector<double> points,
                    const bool do_plot = false);

  template <typename... Args>
  void saveAtPoints(const std::string &filename, bool do_plot, Args... args);

private:
  std::size_t M_Nt;         ///< Number of time steps
  double M_Tfinal;          ///< Final time
  std::vector<dvector> M_T; ///< Solutions in time
  dvector M_xt;             ///< Time discretization
  Model *M_model;           ///< Pointer to a model to solve
};

template <typename... Args>
void SolverTime::saveAtTimes(const std::string &filename, bool do_plot,
                             Args... args) {
  if (M_model->Nx() == 0) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }

  std::filesystem::path filePath(filename);

  // Ensure the directory exists
  std::filesystem::create_directories(filePath.parent_path());

  std::ofstream ofile(filePath, std::ios::out | std::ios::trunc);
  if (!ofile) {
    std::cerr << "Error: can't open " + filename + " !" << std::endl;
    exit(1);
  }

  ofile << "x,";
  ((ofile << args << ","), ...); // Print times to the file

  // Write the solutions for each x at different times
  for (std::size_t i = 0; i < M_model->Nx(); ++i) {
    ofile << "\n" << M_model->x(i);

    ((ofile << "," << M_T[static_cast<std::size_t>(M_Nt * args / M_Tfinal)][i]),
     ...);
  }
  ofile << std::flush;

  if (do_plot) {
    std::string instruction = "python3 ../ploting/plotTimes.py " + filename;
    if (system(instruction.data())) {
      std::cerr << "Error: can't execute plotTimes.py script" << std::endl;
      exit(1);
    }
  }
}

template <typename... Args>
void SolverTime::saveAtPoints(const std::string &filename, bool do_plot,
                              Args... args) {
  if (M_model->Nx() == 0) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }

  std::filesystem::path filePath(filename);

  // Ensure the directory exists
  std::filesystem::create_directories(filePath.parent_path());

  std::ofstream ofile(filePath, std::ios::out | std::ios::trunc);
  if (!ofile) {
    std::cerr << "Error: can't open " + filename + " !" << std::endl;
    exit(1);
  }

  ofile << "x,";
  ((ofile << args << ","), ...); // Print points to the file

  // Write the solutions for each x at different points
  for (std::size_t t = 0; t < M_Nt; ++t) {
    ofile << "\n" << M_xt[t];

    ((ofile << ","
            << M_T[t][static_cast<std::size_t>(M_model->Nx() * args /
                                               M_model->xend())]),
     ...);
  }
  ofile << std::flush;

  if (do_plot) {
    std::string instruction = "python3 ../ploting/plotPoints.py " + filename;
    if (system(instruction.data())) {
      std::cerr << "Error: can't execute plotPoints.py script" << std::endl;
      exit(1);
    }
  }
}

#endif
