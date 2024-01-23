#include "../include/solverTime.hpp"
#include "../include/tridiag.hpp"
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

SolverTime::SolverTime(std::size_t Nt, double tf)
    : M_Nt(Nt), M_Tfinal(tf), M_T(), M_xt(Nt), M_model(nullptr) {}

SolverTime::SolverTime(const SolverTime &s)
    : M_Nt(s.M_Nt), M_Tfinal(s.M_Tfinal), M_T(s.M_T), M_xt(s.M_xt),
      M_model(s.M_model) {}

SolverTime::~SolverTime() {}

double SolverTime::dt() const { return M_Tfinal / M_Nt; }
void SolverTime::setTfinal(double t) { M_Tfinal = t; }
void SolverTime::setNt(std::size_t nt) {
  M_Nt = nt;
  M_xt.resize(nt + 1);
}

void SolverTime::setModel(Model &m) { M_model = &m; }

void SolverTime::solveStatic() {
  if (M_model == nullptr) {
    std::cerr << "Error: you need to set a model first" << std::endl;
    exit(1);
  }

  M_T.clear();
  dvector tmp(M_model->Nx() + 1);

  M_model->setB(tmp, 0, M_Nt, M_Tfinal);
  M_model->setLU();

  solveLU(tmp, M_model->S(), M_model->b());
  M_T.push_back(tmp);
}

void SolverTime::solve() {
  if (M_model == nullptr) {
    std::cerr << "Error: you need to set a model first" << std::endl;
    exit(1);
  }

  M_T.clear();
  dvector tmp(M_model->Nx() + 1);

  M_model->setDt(dt());
  M_model->setLU();

  M_T.push_back(M_model->u()); // Initial conditions in the model
  for (std::size_t t = 0; t < M_Nt + 1; ++t) {
    M_xt[t] = t * dt();

    M_model->setB(M_T.back(), t, M_Nt, M_Tfinal);

    solveLU(tmp, M_model->S(), M_model->b());
    M_T.push_back(tmp);
  }
}

void SolverTime::solve(double tFinal) {
  M_Tfinal = tFinal;
  solve();
}

void SolverTime::saveStatic(const std::string filename, const bool do_plot) {
  if (M_model == nullptr) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }
  if (M_model->Nx() == 0) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }
  std::filesystem::path filePath(filename);

  // Ensure the directory exists
  std::filesystem::create_directories(filePath.parent_path());

  std::ofstream ofile(filePath, std::ios::out | std::ios::trunc);
  if (ofile) {
    ofile << "x,sol,solEx\n";
    for (std::size_t i = 0; i < M_model->Nx(); ++i) {
      ofile << M_model->x(i) << "," << M_T.back()[i] << ","
            << M_model->solveExact(M_model->x(i)) << "\n";
    }
    ofile << std::flush;
  } else {
    std::cerr << "Error: cant open " + filename + " !" << std::endl;
    exit(1);
  }
  std::string instruction = "python3 ../ploting/plotStatic.py " + filename;

  if (do_plot && system(instruction.data())) {
    std::cerr << "Error: cant execut plotStatic.py script " << std::endl;
    exit(1);
  }
}

void SolverTime::saveAtTimes(const std::string filename, dvector times,
                             const bool do_plot) {
  if (M_model == nullptr || M_model->Nx() == 0) {
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
  for (const double &t : times) {
    ofile << t << ",";
  }

  // Write the solutions for each x at different times
  for (std::size_t i = 0; i < M_model->Nx(); ++i) {
    ofile << "\n" << M_model->x(i);

    for (const double &t : times) {
      // Find the corresponding time step index based on the chosen time 't'
      std::size_t index = static_cast<std::size_t>(M_Nt * t / M_Tfinal);

      if (index >= M_T.size()) {
        std::cerr << "Error: invalid time step index" << std::endl;
        exit(1);
      }
      ofile << "," << M_T[index][i];
    }
  }
  if (do_plot) {
    std::string instruction = "python3 ../ploting/plotTimes.py " + filename;
    if (system(instruction.data())) {
      std::cerr << "Error: can't execute plotTimes.py script" << std::endl;
      exit(1);
    }
  }
}

void SolverTime::saveAtPoints(const std::string filename,
                              const std::vector<double> points,
                              const bool do_plot) {
  if (M_model == nullptr || M_model->Nx() == 0) {
    std::cerr << "Error: no data to save" << std::endl;
    exit(1);
  }
  std::ofstream ofile(filename, std::ios::out | std::ios::trunc);
  if (ofile) {
    ofile << "x,";
    for (const double &p : points) {
      ofile << p << ",";
    }

    for (std::size_t t = 0; t < M_Nt; ++t) {
      ofile << "\n" << M_xt[t];

      for (double point : points) {
        std::size_t index =
            static_cast<std::size_t>(M_model->Nx() * point / M_model->xend());
        ofile << "," << M_T[t][index];
      }
    }
    ofile << std::flush;
  } else {
    std::cerr << "Error: can't open " + filename + " !" << std::endl;
    exit(1);
  }

  std::string instruction = "python3 ../ploting/plotPoints.py " + filename;

  if (do_plot && system(instruction.data())) {
    std::cerr << "Error: can't execute plotPoints.py script " << std::endl;
    exit(1);
  }
}
