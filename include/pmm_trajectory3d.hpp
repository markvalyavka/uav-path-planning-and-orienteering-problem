#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "pmm_trajectory.hpp"

#define CONVERGENCE_PRECISION (1E-3)

class PointMassTrajectory3D {
 public:
  PointMassTrajectory3D(const double from, const double to);

 private:
  PMMTrajectory x;
  PMMTrajectory y;
  PMMTrajectory z;
};