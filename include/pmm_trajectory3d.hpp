#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "quad_state.hpp"
#include "pmm_trajectory.hpp"

#define CONVERGENCE_PRECISION (1E-3)

namespace agi {

class PointMassTrajectory3D {
 public:
  PointMassTrajectory3D();
  /*
  version of per axis acceleration limits
  */
  PointMassTrajectory3D(QuadState &from, const QuadState &to,
                        const Vector<3> max_acc,
                        const bool equalize_time = true,
                        const bool calc_gradient = false);
  PointMassTrajectory3D(QuadState &from, const QuadState &to,
                        const Vector<3> max_acc1, const Vector<3> max_acc2,
                        const bool equalize_time = true,
                        const bool calc_gradient = false);
  /*
  version of norm acceleration limits, using GD to distribute acc
  */
  PointMassTrajectory3D(QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm,
                        const bool equalize_time = true,
                        const bool calc_gradient = false);
  /*
  version that limit the thrust by norm using iterative scaling
  */
  PointMassTrajectory3D(QuadState &from, const QuadState &to,
                        const Scalar max_acc_norm, const int max_iter,
                        const bool calc_gradient = false);


  bool exists() const { return x_.exists_ && y_.exists_ && z_.exists_; };
  Scalar time() const { return std::max({x_.time(), y_.time(), z_.time()}); };
  Scalar time_min() const {
    return std::min({x_.time(), y_.time(), z_.time()});
  }
  QuadState state_in_time(const Scalar time_in_tr) const;
  QuadState get_start_state() const;
  QuadState get_end_state() const;
  Vector<3> start_acc() const;
  Vector<3> end_acc() const;
  Vector<3> dt_da_() const {
    return Vector<3>(x_.dt_da_, y_.dt_da_, z_.dt_da_);
  };
  Vector<3> dt_dvs_() const {
    return Vector<3>(x_.dt_dvs_, y_.dt_dvs_, z_.dt_dvs_);
  };
  Vector<3> dt_dve_() const {
    return Vector<3>(x_.dt_dve_, y_.dt_dve_, z_.dt_dve_);
  };

  Vector<3> max_end_velocity_abs() const {
    return Vector<3>(x_.max_end_velocity_abs(), y_.max_end_velocity_abs(),
                     z_.max_end_velocity_abs());
  };

  Scalar get_length_between_times(const Scalar tfrom, const Scalar tto) const;
  Scalar get_length_const_a(const Scalar tfrom, const Scalar tto,
                            const Vector<3> p, const Vector<3> v,
                            const Vector<3> a) const;
  void set_axis_trajectory(const int i, const PMMTrajectory tr);
  PMMTrajectory &get_axis_trajectory(const int i);
  Scalar get_axis_switch_time(const int i) const;
  Scalar inp_from_v_norm;
  Scalar inp_from_v_angle;
  Scalar inp_to_v_norm;
  Scalar inp_to_v_angle;

  PMMTrajectory x_;
  PMMTrajectory y_;
  PMMTrajectory z_;

 private:
  void reproject_to_sphere(Vector<3> &new_thrust_acc, std::vector<bool> &fixed,
                           const Vector<3> &req_thrust_acc_min,
                           const Vector<3> &req_thrust_acc_max,
                           const Vector<3> &acc_req, const Vector<3> &t_times,
                           const Scalar &a_max);
  void get_time_avg_min(Scalar &tavg, Scalar &tmin, int &max_time_idx,
                        int &min_time_idx, const Scalar &tmax,
                        const Vector<3> &t_times,
                        const Vector<3> &gradients_scaled);
};

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &f);
}  // namespace agi