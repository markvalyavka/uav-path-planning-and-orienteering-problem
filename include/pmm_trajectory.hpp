#pragma once
#include <cmath>

#include "agilib/math/types.hpp"

namespace agi {

static constexpr Scalar MAX_SCALAR = std::numeric_limits<Scalar>::max();

class PMMTrajectory {
 public:
  PMMTrajectory();
  PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                const Scalar ve, const Scalar a1_in, const Scalar a2_in,
                const int i, const bool keep_acc_sign = false,
                const bool calc_gradient = false,
                const Scalar check_result = false);
  Scalar time() const { return t_.sum(); }
  Vector<3> state_in_time(const Scalar time_in_tr);
  std::tuple<PMMTrajectory, PMMTrajectory> split_in_time(
    const Scalar time_in_tr);
  void save_to_file(std::string filename);

  bool exists_;
  Vector<3> t_;
  Vector<4> p_;
  Vector<3> v_;
  Vector<2> a_;
  Scalar dt_da_;
  int i_;
};
}  // namespace agi