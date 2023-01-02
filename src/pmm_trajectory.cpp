#include "pmm_trajectory.hpp"

#include <cmath>
#include <fstream>
#include <iostream>


namespace agi {

PMMTrajectory::PMMTrajectory()
  : exists_(false),
    t_(Vector<3>(0, 0, 0)),
    p_(Vector<4>(0, 0, 0, 0)),
    v_(Vector<3>(0, 0, 0)),
    a_(Vector<2>(0, 0)),
    i_(0),
    dt_da_(0) {}

PMMTrajectory::PMMTrajectory(const PMMTrajectory &in) { copy_trajectory(in); }


PMMTrajectory::PMMTrajectory(const Scalar ps, const Scalar vs, const Scalar pe,
                             const Scalar ve, const Scalar a1_in,
                             const Scalar a2_in, const int i,
                             const double desired_time,
                             const bool keep_acc_sign, const bool calc_gradient,
                             const bool check_result) {
  i_ = i;
  p_(0) = ps;
  p_(3) = pe;
  v_(0) = vs;
  v_(2) = ve;

  // std::cout << "\t ps " << ps << " vs " << vs << " pe " << pe << " ve " << ve
  //           << " a1_in " << a1_in << " a2_in " << a2_in << std::endl;

  // If the acceleration is too low, you can't move, and thus, the trajectory
  // doesn't exist.
  if (fabs(a1_in) <= PRECISION_PMM_VALUES &&
      fabs(a2_in) <= PRECISION_PMM_VALUES) {
    exists_ = false;
    return;
  }

  const Scalar pow_vs2 = vs * vs;
  const Scalar pow_ve2 = ve * ve;

  // already there
  if (fabs(pe - ps) < PRECISION_PMM_VALUES &&
      fabs(ve - vs) < PRECISION_PMM_VALUES) {
    a_(0) = 0;
    a_(1) = 0;
    p_(0) = p_(1) = p_(2) = p_(3) = ps;
    v_(0) = v_(1) = v_(2) = vs;
    std::cout << "already there" << std::endl;
    exists_ = true;
    return;
  }

  Scalar t1 = MAX_SCALAR;
  Scalar t2 = MAX_SCALAR;
  Scalar dt_da = MAX_SCALAR;
  Scalar dt_dvs = MAX_SCALAR;
  Scalar dt_dve = MAX_SCALAR;

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

    // std::cout << "\tt1_1 " << t1_1 << " t2_1 " << t2_1 << std::endl;
    // std::cout << "\tt1_2 " << t1_2 << " t2_2 " << t2_2 << std::endl;

    if (std::isfinite(t1_1) and std::isfinite(t2_1) and
        t1_1 > -PRECISION_PMM_VALUES and t2_1 > -PRECISION_PMM_VALUES and
        fabs(desired_time - (t1_1 + t2_1)) < fabs(desired_time - (t1 + t2))) {
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

        const Scalar d_t_dvs = (a1 * vs - a2 * vs - tst1) / (a1 * tst1);
        const Scalar d_t_dve = (-(a1 * ve) + a2 * ve + tst2) / (a2 * tst1);
        // std::cout << i << " d_t_dvs " << d_t_dvs << std::endl;
        // std::cout << i << " d_t_dve " << d_t_dve << std::endl;

        dt_dvs = d_t_dvs;
        dt_dve = d_t_dve;

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
    if (std::isfinite(t1_2) and std::isfinite(t2_2) and
        t1_2 > -PRECISION_PMM_VALUES and t2_2 > -PRECISION_PMM_VALUES and
        (t1 == MAX_SCALAR ||
         fabs(desired_time - (t1_2 + t2_2)) < fabs(desired_time - (t1 + t2)))) {
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


        const Scalar d_t_dvs = (-(a1 * vs) + a2 * vs - tst1) / (a1 * tst1);
        const Scalar d_t_dve = (a1 * ve - a2 * ve + tst2) / (a2 * tst1);

        dt_dvs = d_t_dvs;
        dt_dve = d_t_dve;
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
    dt_dvs_ = dt_dvs;
    dt_dve_ = dt_dve;

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

PMMTrajectory::PMMTrajectory(const PMMTrajectory &in, const Scalar total_time,
                             const bool calc_gradient)
  : PMMTrajectory(in) {
  // recalc the trajectory for known time
  if (time() == total_time) {
    return;
  }
  if (total_time < in.time()) {
    std::cout << "this should never happen" << std::endl;
    std::cout << "in.time() " << in.time() << std::endl;
    std::cout << "total_time " << total_time << std::endl;
  }

  const Scalar &ps = in.p_(0);
  const Scalar &pe = in.p_(3);
  const Scalar &vs = in.v_(0);
  const Scalar &ve = in.v_(2);
  const Scalar a1 = in.a_(0);
  const Scalar a2 = in.a_(1);
  const int &i = in.i_;
  const Scalar pow_ve2 = ve * ve;
  const Scalar pow_vs2 = vs * vs;
  const Scalar pow_a1_2 = a1 * a1;
  const Scalar pow_a2_2 = a2 * a2;
  const Scalar pow_tot2 = total_time * total_time;
  // std::cout << "equalize time " << i << std::endl;

  const Scalar part = (2 * a1 * pe - 2 * a2 * pe - 2 * a1 * ps + 2 * a2 * ps -
                       2 * a1 * total_time * ve + 2 * a2 * total_time * vs);
  const Scalar sqrt_part = sqrt(
    part * part - 4 * a1 * a2 * pow_tot2 * (pow_ve2 - 2 * ve * vs + pow_vs2));
  const Scalar scale1 =
    (-2 * a1 * pe + 2 * a2 * pe + 2 * a1 * ps - 2 * a2 * ps +
     2 * a1 * total_time * ve - 2 * a2 * total_time * vs - sqrt_part) /
    (2 * a1 * a2 * pow_tot2);
  const Scalar scale2 =
    (-2 * a1 * pe + 2 * a2 * pe + 2 * a1 * ps - 2 * a2 * ps +
     2 * a1 * total_time * ve - 2 * a2 * total_time * vs + sqrt_part) /
    (2 * a1 * a2 * pow_tot2);


  // tests
  // const Scalar to_pow = (2 * a1 * pe - 2 * a2 * pe - 2 * a1 * ps + 2 * a2 *
  // ps -
  //                        2 * a1 * total_time * ve + 2 * a2 * total_time *
  //                        vs);
  // const Scalar sqrt_part2 =
  //   sqrt(to_pow * to_pow -
  //        4 * a1 * a2 * pow_tot2 * (pow_ve2 - 2 * ve * vs + pow_vs2));
  // const Scalar t1_1 =
  //   (-(a1 * pe) + a2 * pe + a1 * ps - a2 * ps + a1 * total_time * ve -
  //    a2 * total_time * ve + sqrt_part2 / 2.0) /
  //   (a1 * ve - a2 * ve - a1 * vs + a2 * vs);
  // const Scalar t2_1 =
  //   (-(a1 * pe) + a2 * pe + a1 * ps - a2 * ps + a1 * total_time * vs -
  //    a2 * total_time * vs + sqrt_part2 / 2.0) /
  //   (-(a1 * ve) + a2 * ve + a1 * vs - a2 * vs);
  // const Scalar t1_2 =
  //   (-(a1 * pe) + a2 * pe + a1 * ps - a2 * ps + a1 * total_time * ve -
  //    a2 * total_time * ve - sqrt_part2 / 2.0) /
  //   (a1 * ve - a2 * ve - a1 * vs + a2 * vs);
  // const Scalar t2_2 =
  //   (-(a1 * pe) + a2 * pe + a1 * ps - a2 * ps + a1 * total_time * vs -
  //    a2 * total_time * vs - sqrt_part2 / 2.0) /
  //   (-(a1 * ve) + a2 * ve + a1 * vs - a2 * vs);
  // const Scalar scake1_2 =
  //   (-2 * a1 * pe + 2 * a2 * pe + 2 * a1 * ps - 2 * a2 * ps +
  //    2 * a1 * total_time * ve - 2 * a2 * total_time * vs - sqrt_part2) /
  //   (2 * a1 * a2 * pow_tot2);
  // const Scalar scake2_2 =
  //   (-2 * a1 * pe + 2 * a2 * pe + 2 * a1 * ps - 2 * a2 * ps +
  //    2 * a1 * total_time * ve - 2 * a2 * total_time * vs + sqrt_part2) /
  //   (2 * a1 * a2 * pow_tot2);

  // std::cout << "t1_1 " << t1_1 << " t2_1 " << t2_1 << " sum " << (t1_1 +
  // t2_1)
  //           << std::endl;
  // std::cout << "t1_2 " << t1_2 << " t2_2 " << t2_2 << " sum " << (t1_2 +
  // t2_2)
  //           << std::endl;
  // std::cout << "scake1_2 " << scake1_2 << std::endl;
  // std::cout << "scake2_2 " << scake2_2 << std::endl;
  // std::cout << "ps " << ps << " vs " << vs << std::endl;
  // std::cout << "pe " << pe << " ve " << ve << std::endl;
  // // test times 1
  // const Scalar p1_1 = ps + t1_1 * vs + 0.5 * scake1_2 * a1 * t1_1 * t1_1;
  // const Scalar v1_1 = vs + scake1_2 * a1 * t1_1;

  // const Scalar pe_1 = p1_1 + t2_1 * v1_1 + 0.5 * scake1_2 * a2 * t2_1 * t2_1;
  // const Scalar ve_1 = v1_1 + scake1_2 * a2 * t2_1;
  // std::cout << "pe_1 " << pe_1 << " ve_1 " << ve_1 << std::endl;
  // std::cout << "scake1_2 * a1 " << scake1_2 * a1 << "scake1_2 * a2 "
  //           << scake1_2 * a2 << std::endl;

  // // test times 2
  // const Scalar p1_2 = ps + t1_2 * vs + 0.5 * scake2_2 * a1 * t1_2 * t1_2;
  // const Scalar v1_2 = vs + scake2_2 * a1 * t1_2;

  // const Scalar pe_2 = p1_2 + t2_2 * v1_2 + 0.5 * scake2_2 * a2 * t2_2 * t2_2;
  // const Scalar ve_2 = v1_2 + scake2_2 * a2 * t2_2;
  // std::cout << "pe_2 " << pe_2 << " ve_2 " << ve_2 << std::endl;
  // std::cout << "scake2_2 * a1 " << scake2_2 * a1 << "scake2_2 * a2 "
  //           << scake2_2 * a2 << std::endl;

  // std::cout << std::endl << std::endl;
  // std::cout << "scale1 " << scale1 << std::endl;

  // std::cout << "in.time() " << in.time() << std::endl;

  exists_ = false;
  // std::cout << "total_time " << total_time << std::endl;
  // if (scale1 < 1.0 && scale1 > -1.0) {
  PMMTrajectory scaled =
    PMMTrajectory(ps, vs, pe, ve, a1 * scale1, a2 * scale1, i, total_time, true,
                  calc_gradient, false);
  // std::cout << "1scaled.time() " << scaled.time() << " exists "
  //           << scaled.exists_ << std::endl;
  if (scaled.exists_ and
      fabs(scaled.time() - total_time) <= PRECISION_PMM_VALUES) {
    copy_trajectory(scaled);
  }
  // }
  // std::cout << "scale2 " << scale2 << std::endl;
  if (!exists_) {  //&& scale2 < 1.0 && scale2 > -1.0) {
    PMMTrajectory scaled =
      PMMTrajectory(ps, vs, pe, ve, a1 * scale2, a2 * scale2, i, total_time,
                    true, calc_gradient, false);
    // std::cout << "2scaled.time() " << scaled.time() << " exists "
    //           << scaled.exists_ << std::endl;
    if (scaled.exists_ and
        fabs(scaled.time() - total_time) <= PRECISION_PMM_VALUES) {
      copy_trajectory(scaled);
    }
  }
  if (!exists_) {
    // nummerical issues solution
    // std::cout << "numerical issues scaling" << std::endl;
    // std::cout << "in " << in << std::endl;
    // std::cout << "total_time " << total_time << std::endl;
    // std::cout << "scale1 " << scale1 << std::endl;
    // std::cout << "scale2 " << scale2 << std::endl;
    // exit(1);
    exists_ = false;
    if (fabs(scale1 - 1.0) < 1e-10 || fabs(scale2 - 1.0) < 1e-10) {
      if (t_(0) > 0) {
        t_(0) += total_time - time();
      } else if (t_(2) > 0) {
        t_(2) += total_time - time();
      } else {
        std::cout << "some error here?" << std::endl;
        exists_ = false;
      }
    }
  }

  if (fabs(a_(0)) > fabs(a1)) {
    // std::cout << "above acc" << std::endl;
    exists_ = false;
  }
}

Vector<3> PMMTrajectory::state_in_time(const Scalar time_in_tr) const {
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
    pos = p_(3);
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

void PMMTrajectory::copy_trajectory(const PMMTrajectory &in) {
  exists_ = in.exists_;
  t_ = in.t_;
  p_ = in.p_;
  v_ = in.v_;
  a_ = in.a_;
  i_ = in.i_;
  dt_da_ = in.dt_da_;
  dt_dvs_ = in.dt_dvs_;
  dt_dve_ = in.dt_dve_;
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

Scalar PMMTrajectory::max_end_velocity_abs() const {
  /*what is the end velocity that is when only one acc part is used*/
  Scalar a = a_(0);
  if (v_(2) > v_(0)) {
    if (a < 0) a = a_(1);  // need positive acc
  } else {
    if (a > 0) a = a_(1);  // need negative acc
  }
  const Scalar vs_pow_2 = v_(0) * v_(0);
  Scalar max_v_a = 1.0 / (2.0 * a);
  Scalar max_v_c = p_(0) - p_(3) - vs_pow_2 / (2.0 * a);
  Scalar max_v_disc_sqrt = sqrt(-4 * max_v_a * max_v_c);
  Scalar max_v_abs = max_v_disc_sqrt / (2.0 * max_v_a);
  return max_v_abs;
}

Scalar PMMTrajectory::minRequiredAcc(const Scalar ps, const Scalar vs,
                                     const Scalar pe, const Scalar ve) {
  // required acc to be able to fulfill the task
  // i.e. accelerate between required velocities and in the same time not
  // overshoot the position derived from
  // t_change_vel = (ve-vs)/amax
  // and pe-ps = vs*t_change_vel + 0.5*amax*t_change_vel**2
  // from that what is the amax....
  if (fabs(pe - ps) < PRECISION_TRANS3D) {
    return MIN_ACC_REQ;
  } else {
    const Scalar pow_ve2 = ve * ve;
    const Scalar pow_vs2 = vs * vs;
    if (fabs(pow_ve2 - pow_vs2) < PRECISION_TRANS3D) {
      // does not need any acc basically
      return std::copysign(INITIAL_MIN_ACC, pe - ps);
    } else {
      const Scalar a_min_req = (-pow_ve2 + pow_vs2) / (2 * (-pe + ps));
      return a_min_req;
    }
  }
}

std::ostream &operator<<(std::ostream &o, const PMMTrajectory &f) {
  o << "maxacc: "
    << "t_tot:" << (f.time()) << ";t1:" << f.t_(0) << ";t2:" << f.t_(1)
    << ";t3:" << f.t_(2) << ";exists:" << f.exists_ << ";a1:" << f.a_(0)
    << ";a2:" << f.a_(1) << ";i:" << f.i_ << "\n";
  o << " \t: p0:" << f.p_(0) << ";p1:" << f.p_(1) << ";p2:" << f.p_(2)
    << ";p3:" << f.p_(3) << "\n";
  o << " \t: v0:" << f.v_(0) << ";v1:" << f.v_(1) << ";v3:" << f.v_(2) << "\n";
  o << " \t: dt/dvs:" << f.dt_dvs_ << ";dt/dve:" << f.dt_dve_ << "\n";
  return o;
}

}  // namespace agi