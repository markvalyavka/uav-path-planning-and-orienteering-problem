#include "pmm_trajectory.hpp"

#include <gtest/gtest.h>

#include "agilib/math/types.hpp"

using namespace agi;

Scalar randScalar(Scalar min, Scalar max) {
  Scalar r = (Scalar)rand() / RAND_MAX;
  return min + r * (max - min);
}

/// test random start goals
TEST(SingleAxisPmmTest, RandomStartGoal) {
  static constexpr int N = 1000;

  for (int i = 0; i < N; ++i) {
    // from-to vel in 0-20 m/s
    // pos -50m - 50m
    // a1,a2 in 0--30 / -30--0
    const Scalar vel_from = randScalar(-20, 20);
    const Scalar vel_to = randScalar(-20, 20);
    const Scalar pos_from = randScalar(-50, 50);
    const Scalar pos_to = randScalar(-50, 50);
    const Scalar a1 = randScalar(0, 30);
    const Scalar a2 = randScalar(-30, 0);
    PMMTrajectory tr(pos_from, vel_from, pos_to, vel_to, a1, a2, 0, false, true,
                     true);
    EXPECT_EQ(tr.exists_, true);
  }
}