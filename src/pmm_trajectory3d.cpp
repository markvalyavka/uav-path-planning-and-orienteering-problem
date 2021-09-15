#include "pmm_trajectory3d.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace agi {

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             Scalar max_acc,
                                             const bool equalize_time) {
  x_ =
    PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc, -max_acc, 0);
  y_ =
    PMMTrajectory(from.p(1), from.v(1), to.p(1), to.v(1), max_acc, -max_acc, 1);
  z_ =
    PMMTrajectory(from.p(2), from.v(2), to.p(2), to.v(2), max_acc, -max_acc, 2);

  if (equalize_time) {
    const Scalar tr_time = time();
    for (size_t i = 0; i < 3; i++) {
      if (get_axis_trajectory(i).time() != tr_time) {
        PMMTrajectory scaled = PMMTrajectory(get_axis_trajectory(i), tr_time);
        set_axis_trajectory(i, scaled);
      }
    }
  }
}

QuadState PointMassTrajectory3D::state_in_time(const Scalar time_in_tr) const {
  QuadState ds;
  ds.setZero();
  Vector<3> x = x_.state_in_time(time_in_tr);
  Vector<3> y = y_.state_in_time(time_in_tr);
  Vector<3> z = z_.state_in_time(time_in_tr);

  ds.p(0) = x(0);
  ds.p(1) = y(0);
  ds.p(2) = z(0);

  ds.v(0) = x(1);
  ds.v(1) = y(1);
  ds.v(2) = z(1);

  ds.a(0) = x(2);
  ds.a(1) = y(2);
  ds.a(2) = z(2);

  return ds;
}

QuadState PointMassTrajectory3D::get_start_state() const {
  return state_in_time(0);
}

QuadState PointMassTrajectory3D::get_end_state() const {
  return state_in_time(time());
}

Scalar PointMassTrajectory3D::get_length_between_times(const Scalar tfrom,
                                                       const Scalar tto) const {
  // time switches in t1 for x,y,z
  std::vector<Scalar> switch_times = {tfrom, tto};
  for (size_t i = 0; i < 3; i++) {
    const Scalar switch_time = get_axis_switch_time(i);
    if (switch_time > tfrom and switch_time < tto) {
      switch_times.push_back(switch_time);
    }
  }

  // sort the times
  std::sort(switch_times.begin(), switch_times.end());

  Scalar ds = 0;
  for (size_t i = 1; i < switch_times.size(); i++) {
    const Scalar tfrom_part = switch_times[i - 1];
    const Scalar tto_part = switch_times[i];
    QuadState qs = state_in_time(tfrom_part);

    ds += get_length_const_a(0, tto_part - tfrom_part, qs.p, qs.v, qs.a);
  }
  return ds;
}

Scalar PointMassTrajectory3D::get_length_const_a(const Scalar tfrom,
                                                 const Scalar tto,
                                                 const Vector<3> p,
                                                 const Vector<3> v,
                                                 const Vector<3> a) const {
  // const double &sx = p(0);
  // const double &sy = p(1);
  // const double &sz = p(2);

  const double ax_pow2 = a(0) * a(0);
  const double ay_pow2 = a(1) * a(1);
  const double az_pow2 = a(2) * a(2);
  const double vx_pow2 = v(0) * v(0);
  const double vy_pow2 = v(1) * v(1);
  const double vz_pow2 = v(2) * v(2);

  double t = tto;
  double t_pow2 = t * t;

  double tmp = sqrt(ax_pow2 * t_pow2 + ay_pow2 * t_pow2 + az_pow2 * t_pow2 +
                    2 * a(0) * t * v(0) + vx_pow2 + 2 * a(1) * t * v(1) +
                    vy_pow2 + 2 * a(2) * t * v(2) + vz_pow2);
  double logpart =
    log(ax_pow2 * t + ay_pow2 * t + az_pow2 * t + a(0) * v(0) + a(1) * v(1) +
        a(2) * v(2) + sqrt(ax_pow2 + ay_pow2 + az_pow2) * tmp);
  if (v(0) == 0 && v(1) == 0 && v(2) == 0 && t == 0) {
    logpart = 0;
  }
  const double ds_to =
    (0.5) *
    (tmp * (t + (a(0) * v(0) + a(1) * v(1) + a(2) * v(2)) /
                  (ax_pow2 + ay_pow2 + az_pow2)) +
     (1.0 / pow(ax_pow2 + ay_pow2 + az_pow2, 1.5)) *
       ((az_pow2 * (vx_pow2 + vy_pow2) - 2.0 * a(0) * a(2) * v(0) * v(2) -
         2.0 * a(1) * v(1) * (a(0) * v(0) + a(2) * v(2)) +
         ay_pow2 * (vx_pow2 + vz_pow2) + ax_pow2 * (vy_pow2 + vz_pow2)) *
        logpart));


  t = tfrom;
  t_pow2 = t * t;

  tmp = sqrt(ax_pow2 * t_pow2 + ay_pow2 * t_pow2 + az_pow2 * t_pow2 +
             2 * a(0) * t * v(0) + vx_pow2 + 2 * a(1) * t * v(1) + vy_pow2 +
             2 * a(2) * t * v(2) + vz_pow2);
  logpart =
    log(ax_pow2 * t + ay_pow2 * t + az_pow2 * t + a(0) * v(0) + a(1) * v(1) +
        a(2) * v(2) + sqrt(ax_pow2 + ay_pow2 + az_pow2) * tmp);
  if (v(0) == 0 && v(1) == 0 && v(2) == 0 && t == 0) {
    logpart = 0;
  }

  const double ds_from =
    (0.5) *
    (tmp * (t + (a(0) * v(0) + a(1) * v(1) + a(2) * v(2)) /
                  (ax_pow2 + ay_pow2 + az_pow2)) +
     (1.0 / pow(ax_pow2 + ay_pow2 + az_pow2, 1.5)) *
       ((az_pow2 * (vx_pow2 + vy_pow2) - 2.0 * a(0) * a(2) * v(0) * v(2) -
         2.0 * a(1) * v(1) * (a(0) * v(0) + a(2) * v(2)) +
         ay_pow2 * (vx_pow2 + vz_pow2) + ax_pow2 * (vy_pow2 + vz_pow2)) *
        logpart));

  const double ds = ds_to - ds_from;

  if (!std::isfinite(ds) || ds < 0) {
    std::cout << "non finite ds or bellow zero" << std::endl;
    std::cout << "p " << p.transpose() << std::endl;
    std::cout << "v " << v.transpose() << std::endl;
    std::cout << "a " << a.transpose() << std::endl;
    return INF;
  }
  return ds;
}

void PointMassTrajectory3D::set_axis_trajectory(const int i,
                                                const PMMTrajectory tr) {
  switch (i) {
    case 0:
      x_ = tr;
      break;
    case 1:
      y_ = tr;
      break;
    case 2:
      z_ = tr;
      break;
    default:
      exit(1);
  }
}
PMMTrajectory &PointMassTrajectory3D::get_axis_trajectory(const int i) {
  switch (i) {
    case 0:
      return x_;
    case 1:
      return y_;
    case 2:
      return z_;
    default:
      exit(1);
  }
}

Scalar PointMassTrajectory3D::get_axis_switch_time(const int i) const {
  switch (i) {
    case 0:
      return x_.t_(0);
    case 1:
      return y_.t_(0);
    case 2:
      return z_.t_(0);
    default:
      exit(1);
  }
}

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &t) {
  o << "pmm3d: t:" << t.time() << ";exists:" << t.exists();
  o << "\n\tx: " << t.x_;
  o << "\n\ty: " << t.y_;
  o << "\n\tz: " << t.z_;
  return o;
}

}  // namespace agi