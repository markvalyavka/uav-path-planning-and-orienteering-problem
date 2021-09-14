#include "pmm_trajectory.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#define PRECISION_PMM_VALUES (1.0e-8)

namespace agi {

PMMTrajectory::PMMTrajectory()
  : exists_(false),
    t_(Vector<3>(0, 0, 0)),
    p_(Vector<4>(0, 0, 0, 0)),
    v_(Vector<3>(0, 0, 0)),
    a_(Vector<2>(0, 0)),
    i_(0),
    dt_da_(0) {}

PMMTrajectory::PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                             const Scalar ve, const Scalar a1_in,
                             const Scalar a2_in, const int i,
                             const bool keep_acc_sign, const bool calc_gradient,
                             const Scalar check_result) {
  i_ = i;
  p_(0) = ps;
  p_(3) = pe;
  v_(0) = vs;
  v_(2) = ve;

  // can not move without acc
  if (fabs(a1_in) <= PRECISION_PMM_VALUES &&
      fabs(a2_in) <= PRECISION_PMM_VALUES) {
    exists_ = false;
    return;
  }

  const Scalar pow_ve2 = ve * ve;
  const Scalar pow_vs2 = vs * vs;

  // already there
  if (fabs(pe - ps) < PRECISION_PMM_VALUES &&
      fabs(ve - vs) < PRECISION_PMM_VALUES) {
    a_(0) = a1_in;
    a_(1) = a2_in;
    p_(0) = p_(1) = p_(2) = p_(3) = ps;
    v_(0) = v_(1) = v_(3) = vs;
    exists_ = true;
    return;
  }

  Scalar t1 = MAX_SCALAR;
  Scalar t2 = MAX_SCALAR;
  Scalar dt_da = MAX_SCALAR;

  Scalar used_acc1 = a1_in;
  Scalar used_acc2 = a2_in;
  std::vector<Vector<2>> test_acc_vec = {Vector<2>(a1_in, a2_in),
                                         Vector<2>(a2_in, a1_in)};
  if (keep_acc_sign) {
    test_acc_vec = {Vector<2>(a1_in, a2_in)};
  }
  for (Vector<2> a : test_acc_vec) {
    const Scalar a1 = a(0);
    const Scalar a2 = a(1);
    const Scalar pow_a2_2 = a2 * a2;
    const Scalar pow_a1_2 = a1 * a1;
    const Scalar pow_a1_3 = pow_a1_2 * a1;

    const Scalar tst1 = sqrt(
      (-a1 + a2) * (2 * a1 * a2 * (pe - ps) - a1 * pow_ve2 + a2 * pow_vs2));
    const Scalar tst2 = sqrt(
      (a1 - a2) * (a1 * (-2 * a2 * pe + 2 * a2 * ps + pow_ve2) - a2 * pow_vs2));

    // case 1
    const Scalar t1_1 = (-(a1 * vs) + a2 * vs + tst1) / (a1 * (a1 - a2));
    const Scalar t2_1 = -((-(a1 * ve) + a2 * ve + tst2) / ((a1 - a2) * a2));

    // case 2
    const Scalar t1_2 = -((a1 * vs - a2 * vs + tst1) / (a1 * (a1 - a2)));
    const Scalar t2_2 = (a1 * ve - a2 * ve + tst2) / ((a1 - a2) * a2);

    if (std::isfinite(t1_1) and std::isfinite(t2_1) and
        t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES and
        t1_1 + t2_1 < t1 + t2) {
      t1 = std::max(t1_1, 0.0);
      t2 = std::max(t2_1, 0.0);

      if (calc_gradient) {
        // dt/da == gradient we can use to optimize time
        const Scalar d_t_da1 = (a1 * (2 * a2 * (pe - ps) - pow_ve2 - pow_vs2) +
                                2 * vs * (a2 * vs + tst1)) /
                               (2 * pow_a1_2 * tst1);
        const Scalar d_t_da2 = (2 * a1 * (a2 * (-pe + ps) + pow_ve2) -
                                a2 * (pow_ve2 + pow_vs2) - 2 * ve * tst2) /
                               (2 * pow_a2_2 * tst1);

        if (i < 2) {
          dt_da = std::copysign(d_t_da1, -a1) +
                  std::copysign(d_t_da2, -a1);  // gradient with respect to a1
        } else {
          dt_da = std::copysign(d_t_da1, -1) + std::copysign(d_t_da2, -1);
          // gradient with respect to both a1 and a2 as for z axis those are
          // different
          // a1 always positive
        }
      }
      used_acc1 = a1;
      used_acc2 = a2;
    } else if (std::isfinite(t1_2) and std::isfinite(t2_2) and
               t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES and
               t1_2 + t2_2 < t1 + t2) {
      t1 = std::max(t1_2, 0.0);
      t2 = std::max(t2_2, 0.0);

      if (calc_gradient) {
        // dt/da == gradient we can use to optimize time
        const Scalar d_t_da1 =
          (a1 * (-2 * a2 * pe + 2 * a2 * ps + pow_ve2 + pow_vs2) +
           2 * vs * (-(a2 * vs) + tst1)) /
          (2 * pow_a1_2 * tst1);
        const Scalar d_t_da2 = (2 * a1 * a2 * (pe - ps) - 2 * a1 * pow_ve2 +
                                a2 * (pow_ve2 + pow_vs2) - 2 * ve * tst2) /
                               (2 * pow_a2_2 * tst1);

        if (i < 2) {
          dt_da = std::copysign(d_t_da1, -a1) +
                  std::copysign(d_t_da2, -a1);  // gradient with respect to a1
        } else {
          dt_da = std::copysign(d_t_da1, -1) + std::copysign(d_t_da2, -1);
          // gradient with respect to both a1 and a2 as for z axis those are
          // different
          // a1 always positive
        }
      }
      used_acc1 = a1;
      used_acc2 = a2;
    }
  }

  if (t1 >= 0 and t1 != MAX_SCALAR) {
    exists_ = true;
    t_ = Vector<3>(t1, 0, t2);
    p_(1) = ps + t1 * vs + 0.5 * used_acc1 * t1 * t1;
    p_(2) = p_(1);
    v_(1) = vs + used_acc1 * t1;
    a_(0) = used_acc1;
    a_(1) = used_acc2;
    dt_da_ = std::isfinite(dt_da) ? dt_da : 0;

    if (check_result) {
      const Scalar ve_tst = v_(1) + a_(1) * t2;
      const Scalar pe_tst = p_(2) + t2 * v_(1) + 0.5 * a_(1) * t2 * t2;
      if (fabs(ve_tst - ve) > PRECISION_PMM_VALUES ||
          fabs(pe_tst - pe) > PRECISION_PMM_VALUES) {
        std::cout << "wrong ve or pe oddi two acc" << std::endl;
        std::cout << "ve_tst " << ve_tst << " ve " << ve << std::endl;
        std::cout << "pe_tst " << pe_tst << " pe " << pe << std::endl;
        std::cout << "t1 " << t1 << " t2 " << t2 << " a1 " << used_acc1
                  << " a2 " << used_acc2 << std::endl;
        exists_ = false;
      }
    }

  } else {
    exists_ = false;
  }
}

Vector<3> PMMTrajectory::state_in_time(const Scalar time_in_tr) {
  Scalar pos, vel, acc;

  if (time_in_tr < t_(0)) {  // const acc part with a1
    pos = p_(0) + v_(0) * time_in_tr + 0.5 * a_(0) * time_in_tr * time_in_tr;
    vel = v_(0) + a_(0) * time_in_tr;
    acc = a_(0);
  } else if (time_in_tr < t_(0) + t_(1)) {  // const vel part with a=0
    const Scalar time_part = (time_in_tr - t_(0));
    pos = p_(1) + v_(1) * time_part;
    vel = v_(1);
    acc = 0.0;
  } else if (time_in_tr < time()) {  // const vel part with a2
    const Scalar time_part = (time_in_tr - t_(0) - t_(1));
    pos = p_(2) + v_(1) * time_part + 0.5 * a_(1) * time_part * time_part;
    vel = v_(1) + a_(1) * time_part;
    acc = a_(1);
  } else {  // return the last state
    pos = p_(4);
    vel = v_(2);
    acc = a_(1);
  }

  if (t_(0) == 0 && t_(1) == 0 && t_(2) > 0) {
    acc = a_(1);
  } else if (t_(0) > 0 and t_(1) == 0 and t_(2) == 0) {
    acc = a_(0);
  }

  return Vector<3>(pos, vel, acc);
}

std::tuple<PMMTrajectory, PMMTrajectory> PMMTrajectory::split_in_time(
  const Scalar time_in_tr) {
  Vector<3> state = state_in_time(time_in_tr);

  PMMTrajectory bef;
  bef.exists_ = true;
  bef.i_ = i_;
  bef.a_ = a_;
  PMMTrajectory aft;
  aft.exists_ = true;
  aft.i_ = i_;
  aft.a_ = a_;

  if (time_in_tr <= t_(0)) {
    bef.t_(0) = time_in_tr;
    bef.t_(1) = bef.t_(2) = .0;
    bef.t_(2) = .0;
    bef.p_(0) = p_(0);
    bef.p_(1) = bef.p_(2) = bef.p_(3) = state(0);
    bef.v_(0) = v_(0);
    bef.v_(1) = bef.v_(2) = state(1);

    aft.t_(0) = t_(0) - time_in_tr;
    aft.t_(1) = t_(1);
    aft.t_(2) = t_(2);
    aft.p_(0) = state(0);
    aft.p_(1) = p_(1);
    aft.p_(2) = p_(2);
    aft.p_(3) = p_(3);
    aft.v_(0) = state(1);
    aft.v_(1) = v_(1);
    aft.v_(2) = v_(2);
  } else if (time_in_tr <= t_(0) + t_(1)) {
    bef.t_(0) = t_(0);
    bef.t_(1) = time_in_tr - t_(0);
    bef.t_(2) = .0;
    bef.p_(0) = p_(0);
    bef.p_(1) = p_(1);
    bef.p_(2) = bef.p_(3) = state(0);
    bef.v_(0) = v_(0);
    bef.v_(1) = v_(1);
    bef.v_(2) = state(1);

    aft.t_(0) = .0;
    aft.t_(1) = (t_(0) + t_(1)) - time_in_tr;
    aft.t_(2) = t_(2);
    aft.p_(0) = aft.p_(1) = state(0);
    aft.p_(2) = p_(2);
    aft.p_(3) = p_(3);
    aft.v_(0) = aft.v_(1) = state(1);
    aft.v_(2) = v_(2);
  } else if (time_in_tr <= t_(0) + t_(1) + t_(2)) {
    bef.t_ = t_;
    bef.t_(2) = time_in_tr - (t_(0) + t_(1));
    bef.p_ = p_;
    bef.p_(3) = state(0);
    bef.v_ = v_;
    bef.v_(2) = state(1);

    aft.t_(0) = aft.t_(1) = .0;
    aft.t_(2) = (t_(0) + t_(1) + t_(2)) - time_in_tr;
    aft.p_(0) = aft.p_(1) = aft.p_(2) = state(0);
    aft.p_(3) = p_(3);
    aft.v_(0) = aft.v_(1) = state(1);
    aft.v_(2) = v_(2);
  } else {
    bef.t_ = t_;
    bef.p_ = p_;
    bef.v_ = v_;

    aft.t_(0) = aft.t_(1) = aft.t_(2) = .0;
    aft.p_(0) = aft.p_(1) = aft.p_(3) = state(0);
    aft.v_(0) = aft.v_(1) = aft.v_(2) = state(1);
  }
  return {bef, aft};
}


void PMMTrajectory::save_to_file(std::string filename) {
  std::ofstream myfile;
  myfile.open(filename.c_str());
  if (myfile.is_open()) {
    myfile << "a1,a2,t1,t2,t3,p0,p1,p2,p3,v0,v1,"
              "v3,exists,dt_da,axis"
           << std::endl;
    myfile << a_(0) << "," << a_(1) << "," << t_(0) << "," << t_(1) << ","
           << t_(2) << "," << p_(0) << "," << p_(1) << "," << p_(2) << ","
           << p_(3) << "," << v_(0) << "," << v_(1) << "," << v_(2) << ","
           << exists_ << "," << dt_da_ << "," << i_;
    myfile.close();
  }
}

}  // namespace agi