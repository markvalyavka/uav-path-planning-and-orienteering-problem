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
                        double max_acc);
  bool exists(){return x_.exists_ && y_.exists_ && z_.exists_};
  Scalar time() { return std::max({x_.time(), y_.time(), z_.time()}); };

 private:
  PMMTrajectory x_;
  PMMTrajectory y_;
  PMMTrajectory z_;
};

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &f)
}  // namespace agi