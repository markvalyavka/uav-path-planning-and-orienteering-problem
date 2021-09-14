#pragma once
#include <cmath>
#include <eigen3/Eigen/Eigen>

using Scalar = double;
static constexpr Scalar INF = std::numeric_limits<Scalar>::infinity();
static constexpr Scalar MAX_SCALAR = std::numeric_limits<Scalar>::max();

// Define `Dynamic` matrix size.
static constexpr int Dynamic = Eigen::Dynamic;

// Using shorthand for `Matrix<rows, cols>` with scalar type.
template<int rows = Dynamic, int cols = rows>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

// Using shorthand for `Vector<rows>` with scalar type.
template<int rows = Dynamic>
using Vector = Matrix<rows, 1>;

class PMMTrajectory {
 public:
  PMMTrajectory();
  PMMTrajectory(const double ps, const double vs, const double pe,
                const double ve, const double a1_in, const double a2_in,
                const int i, const bool keep_acc_sign, const bool calc_gradient,
                const double check_result);
  double time() const { return t_.sum(); }
  Vector<3> state_in_time(const Scalar time_in_tr);
  std::tuple<PMMTrajectory, PMMTrajectory> split_in_time(
    const Scalar time_in_tr);
  void save_to_file(std::string filename);

  bool exists_;
  Vector<3> t_;
  Vector<4> p_;
  Vector<3> v_;
  Vector<2> a_;
  double dt_da_;
  int i_;
};