#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "agilib/types/quadrotor.hpp"
#include "pmm_trajectory.hpp"

#define CONVERGENCE_PRECISION (1E-3)

namespace agi {

class PointMassTrajectory3D {
 public:
  PointMassTrajectory3D(const QuadState &from, const QuadState &to,
                        Scalar max_acc, const bool equalize_time = true);
  bool exists() const { return x_.exists_ && y_.exists_ && z_.exists_; };
  Scalar time() const { return std::max({x_.time(), y_.time(), z_.time()}); };
  Scalar time_min() const {
    return std::min({x_.time(), y_.time(), z_.time()});
  }
  QuadState state_in_time(const Scalar time_in_tr) const;
  QuadState get_start_state() const;
  QuadState get_end_state() const;

  Scalar get_length_between_times(const Scalar tfrom, const Scalar tto) const;
  Scalar get_length_const_a(const Scalar tfrom, const Scalar tto,
                            const Vector<3> p, const Vector<3> v,
                            const Vector<3> a) const;
  void set_axis_trajectory(const int i, const PMMTrajectory tr);
  PMMTrajectory &get_axis_trajectory(const int i);
  Scalar get_axis_switch_time(const int i) const;

  PMMTrajectory x_;
  PMMTrajectory y_;
  PMMTrajectory z_;
};

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &f);
}  // namespace agi