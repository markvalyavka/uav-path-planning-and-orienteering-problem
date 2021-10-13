#pragma once
#include <cmath>

#include "agilib/math/types.hpp"

#define PRECISION_PMM_VALUES (1.0e-8)
#define PRECISION_TRANS3D (1E-4)
#define MIN_ACC_REQ (0.01)
#define INITIAL_MIN_ACC (0.5)

namespace agi {

static constexpr Scalar MAX_SCALAR = std::numeric_limits<Scalar>::max();

class PMMTrajectory {
 public:
  PMMTrajectory();
  PMMTrajectory(const PMMTrajectory &in);
  /*
  this solves the equations
  v1 == vs + t1*(a1);
  ve == v1 + t2*(a2);
  p1 == ps + t1*vs + 1/2*t1^2*(a1);
  pe == p1 + t2*v1 + 1/2*t2^2*(a2);
  */
  PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                const Scalar ve, const Scalar a1_in, const Scalar a2_in,
                const int i, const bool keep_acc_sign = false,
                const bool calc_gradient = false,
                const Scalar check_result = false);
  /*
  one that scales times to given one
  */
  PMMTrajectory(const PMMTrajectory &in, const Scalar total_time,
                const double calc_gradient = false);

  static Scalar minRequiredAcc(const Scalar ps, const Scalar vs,
                               const Scalar pe, const Scalar ve);

  Scalar time() const { return t_.sum(); }
  Vector<3> state_in_time(const Scalar time_in_tr) const;
  std::tuple<PMMTrajectory, PMMTrajectory> split_in_time(
    const Scalar time_in_tr);
  void copy_trajectory(const PMMTrajectory &in);
  void save_to_file(std::string filename);

  bool exists_;
  Vector<3> t_;
  Vector<4> p_;
  Vector<3> v_;
  Vector<2> a_;
  Scalar dt_da_;
  Scalar dt_dvs_;
  Scalar dt_dve_;
  int i_;
};

std::ostream &operator<<(std::ostream &o, const PMMTrajectory &f);
}  // namespace agi