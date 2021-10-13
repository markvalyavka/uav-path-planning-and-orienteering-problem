#include "pmm_trajectory3d.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

#include "agilib/math/gravity.hpp"

namespace agi {
PointMassTrajectory3D::PointMassTrajectory3D() {}

/*
basic version with simetric acc limits in axis
*/
PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Vector<3> max_acc,
                                             const bool equalize_time,
                                             const bool calc_gradient)
  : PointMassTrajectory3D(from, to, max_acc, -max_acc, equalize_time,
                          calc_gradient) {}


//                                              {
//   x_ = PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc(0),
//                      -max_acc(0), 0);
//   y_ = PMMTrajectory(from.p(1), from.v(1), to.p(1), to.v(1), max_acc(1),
//                      -max_acc(1), 1);
//   z_ = PMMTrajectory(from.p(2), from.v(2), to.p(2), to.v(2), max_acc(2),
//                      -max_acc(2), 2);


//   if (equalize_time) {
//     const Scalar tr_time = time();
//     for (size_t i = 0; i < 3; i++) {
//       if (get_axis_trajectory(i).time() != tr_time) {
//         PMMTrajectory scaled = PMMTrajectory(get_axis_trajectory(i),
//         tr_time); set_axis_trajectory(i, scaled);
//       }
//     }
//   }
// }

/*
version that converges to thust limit by iterative increasing scaled (to mach
time) to the acc norm
this version scales time by default
*/
PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Scalar max_acc_norm,
                                             const int max_iter,
                                             const bool calc_gradient) {
  // B = A + GVEC , |A|=max_acc_norm
  // initially B equal per axis with b_x=b_y=b_z -> |B-GVEC|^2 = |T|^2
  // -> 3*bs_x^2 + 2*g*a_x + g^2 - |T|^2 = 0 --> roots are the possible acc

  const Scalar precision_acc_limit = 0.1;

  const Scalar max_acc_norm_pow2 = max_acc_norm * max_acc_norm;
  const Scalar a_equal_acc = 3;
  const Scalar b_equal_acc = 2 * G;
  const Scalar c_equal_acc = G * G - max_acc_norm_pow2;
  const Scalar equal_acc_1 =
    (-b_equal_acc +
     sqrt(b_equal_acc * b_equal_acc - 4 * a_equal_acc * c_equal_acc)) /
    (2 * a_equal_acc);
  // const Scalar equal_acc_2 =
  //   (-b_equal_acc -
  //    sqrt(b_equal_acc * b_equal_acc - 4 * a_equal_acc * c_equal_acc)) /
  //   (2 * a_equal_acc);
  // std::cout << "equal_acc_1 " << equal_acc_1 << std::endl;
  // std::cout << "equal_acc_2 " << equal_acc_2 << std::endl;
  // Vector<3> T1_1 = Vector<3>::Constant(equal_acc_1) - GVEC;
  // Vector<3> T1_2 = Vector<3>::Constant(-equal_acc_1) - GVEC;
  // Vector<3> T2_1 = Vector<3>::Constant(equal_acc_2) - GVEC;
  // Vector<3> T2_2 = Vector<3>::Constant(-equal_acc_2) - GVEC;
  // std::cout << "T1_1 " << T1_1.transpose() << " norm " << T1_1.norm()
  //           << std::endl;
  // std::cout << "T1_2 " << T1_2.transpose() << " norm " << T1_2.norm()
  //           << std::endl;
  // std::cout << "T2_1 " << T2_1.transpose() << " norm " << T2_1.norm()
  //           << std::endl;
  // std::cout << "T2_2 " << T2_2.transpose() << " norm " << T2_2.norm()
  //           << std::endl;


  Vector<3> equal_acc;
  if ((Vector<3>::Constant(-equal_acc_1) - GVEC).norm() < max_acc_norm) {
    equal_acc = Vector<3>::Constant(equal_acc_1);
  } else {
    const Scalar equal_acc_2 =
      (-b_equal_acc -
       sqrt(b_equal_acc * b_equal_acc - 4 * a_equal_acc * c_equal_acc)) /
      (2 * a_equal_acc);
    equal_acc = Vector<3>::Constant(equal_acc_2);
    std::cout << "that is happening ever?" << std::endl;
  }


  PointMassTrajectory3D pmm3d(from, to, equal_acc, true);
  // if(!pmm3d.exists())
  Vector<3> start_acc = pmm3d.start_acc();
  Vector<3> end_acc = pmm3d.end_acc();
  Scalar start_thrust = (start_acc - GVEC).norm();
  Scalar end_thrust = (end_acc - GVEC).norm();

  Vector<3> bigger_acc = start_thrust > end_thrust ? start_acc : end_acc;
  Scalar larger_thrus = std::max(start_thrust, end_thrust);
  int iter = 0;

  // std::cout << "start_acc " << start_acc.transpose() << std::endl;
  // std::cout << "end_acc " << end_acc.transpose() << std::endl;
  // std::cout << "start_thrust " << (start_acc - GVEC).transpose() <<
  // std::endl; std::cout << "end_thrust " << (end_acc - GVEC).transpose() <<
  // std::endl;

  while (fabs(larger_thrus - max_acc_norm) > precision_acc_limit and
         iter < max_iter) {
    iter++;
    // B = T + GVEC , |T|=max_acc_norm
    // scale the A parts by same factor k ->
    // k^2*b_x^2 + k^2*b_y^2 + (k*b_z + g)^2 - |T|^2 = 0 -> find k
    // (b_x^2 + b_y^2 + b_z^2) * k^2 + (2*g*b_z)*k


    const Scalar a_dis = bigger_acc.squaredNorm();
    const Scalar b_dis = 2 * bigger_acc(2) * G;
    const Scalar c_dis = G * G - (max_acc_norm * max_acc_norm);
    const Scalar k1 =
      (-b_dis + sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    const Scalar k2 =
      (-b_dis - sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    Vector<3> thrust_acc_new_k1 = k1 * bigger_acc - GVEC;
    Vector<3> thrust_acc_new_k2 = k2 * bigger_acc - GVEC;
    Vector<3> max_acc_new1 = thrust_acc_new_k1 + GVEC;
    Vector<3> max_acc_new2 = -thrust_acc_new_k1 + GVEC;


    pmm3d = PointMassTrajectory3D(from, to, max_acc_new1, max_acc_new2, true,
                                  calc_gradient);
    start_acc = pmm3d.start_acc();
    end_acc = pmm3d.end_acc();
    start_thrust = (start_acc - GVEC).norm();
    end_thrust = (end_acc - GVEC).norm();
    // std::cout << "-----------------------------" << std::endl;
    // std::cout << "start_acc " << start_acc.transpose() << std::endl;
    // std::cout << "end_acc " << end_acc.transpose() << std::endl;
    // std::cout << "start_thrust " << (start_acc - GVEC).transpose() <<
    // std::endl; std::cout << "end_thrust " << (end_acc - GVEC).transpose() <<
    // std::endl;

    if (start_thrust > end_thrust) {
      bigger_acc = start_acc;
      larger_thrus = start_thrust;
    } else {
      bigger_acc = end_acc;
      larger_thrus = end_thrust;
    }

    // std::cout << "acc_vec " << acc_vec.transpose() << std::endl;
    // std::cout << "max_acc_new1 " << max_acc_new1.transpose() << std::endl;
    // std::cout << "max_acc_new2 " << max_acc_new2.transpose() << std::endl;
    // std::cout << "k1 " << k1 << " k1 * acc_vec " << (k1 *
    // acc_vec).transpose()
    //           << " thrust_acc_new_k1 " << thrust_acc_new_k1.transpose()
    //           << " norm " << thrust_acc_new_k1.norm() << std::endl;
    // std::cout << "k2 " << k2 << " k2 * acc_vec " << (k2 *
    // acc_vec).transpose()
    //           << " thrust_acc_new_k2 " << thrust_acc_new_k2.transpose()
    //           << " norm " << thrust_acc_new_k2.norm() << std::endl;
    // std::cout << "pmm3d iter " << iter << std::endl;
    // std::cout << pmm3d << std::endl;
    // std::cout << "pmm start acc" << pmm3d.get_start_state().a.norm()
    //           << " thrust norm " << (pmm3d.get_start_state().a - GVEC).norm()
    //           << std::endl;
    // std::cout << "pmm end acc" << pmm3d.get_end_state().a.norm()
    //           << " thrust norm " << (pmm3d.get_end_state().a - GVEC).norm()
    //           << std::endl;
    // std::cout << "iter " << iter << std::endl;
  }

  x_ = pmm3d.x_;
  y_ = pmm3d.y_;
  z_ = pmm3d.z_;
}

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Vector<3> max_acc1,
                                             const Vector<3> max_acc2,
                                             const bool equalize_time,
                                             const bool calc_gradient) {
  x_ = PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc1(0),
                     max_acc2(0), 0);
  y_ = PMMTrajectory(from.p(1), from.v(1), to.p(1), to.v(1), max_acc1(1),
                     max_acc2(1), 1);
  z_ = PMMTrajectory(from.p(2), from.v(2), to.p(2), to.v(2), max_acc1(2),
                     max_acc2(2), 2);


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

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Scalar max_acc_norm,
                                             const bool equalize_time,
                                             const bool calc_gradient) {
  // std::cout << "gd optimization " << std::endl;
  static const Scalar ALLOWED_DIFF_TIMES_RATIO{0.0001};
  static const Scalar NUM_ITERS{10};

  Vector<3> t_times(0.0, 0.0, 0.0);
  Vector<3> gradients(0.0, 0.0, 0.0);
  Vector<3> acc_req(MIN_ACC_REQ, MIN_ACC_REQ, MIN_ACC_REQ);

  for (int i = 0; i < 3; ++i) {
    const Scalar vs = from.v(i);
    const Scalar ve = to.v(i);
    const Scalar ps = from.p(i);
    const Scalar pe = to.p(i);
    acc_req(i) = PMMTrajectory::minRequiredAcc(ps, vs, pe, ve);
    if (acc_req(i) < 0) {
      acc_req(i) = std::min(acc_req(i), -MIN_ACC_REQ);
    } else {
      acc_req(i) = std::max(acc_req(i), MIN_ACC_REQ);
    }


    PMMTrajectory tr_above(ps, vs, pe, ve, acc_req(i) * 1.1, -acc_req(i) * 1.1,
                           i, false, true, false);
    if (tr_above.exists_) {
      acc_req(i) = std::copysign(acc_req(i), tr_above.a_(0));
      t_times(i) = tr_above.time();
      gradients(i) = tr_above.dt_da_;
    } else {
      std::cout << "non existing min acc should not happen i:" << i
                << std::endl;
      exit(1);
    }
    PMMTrajectory tr_bellow(ps, vs, pe, ve, acc_req(i) * 0.9, -acc_req(i) * 0.9,
                            i, false, true, false);
    if (tr_bellow.exists_) {
      acc_req(i) = std::copysign(MIN_ACC_REQ, acc_req(i));
    } else {
      std::cout << "bellow does not exists i:" << i << std::endl;
    }
  }

  acc_req(2) = fabs(acc_req(2));

  // we need abs of the az + g to be within the limit....
  const Vector<3> req_max_thrust_acc = acc_req.cwiseAbs() + Vector<3>(0, 0, G);

  if (req_max_thrust_acc.norm() > max_acc_norm) {
    x_.exists_ = false;
    y_.exists_ = false;
    z_.exists_ = false;
    std::cout << "req_max_thrust_acc above limit " << req_max_thrust_acc.norm()
              << " a_max " << max_acc_norm << std::endl;
    return;
  }

  // add G to min thrust acc
  Vector<3> req_thrust_acc_min;
  Vector<3> req_thrust_acc_max;
  for (size_t i = 0; i < 3; i++) {
    if (i == 2) {
      // z
      if (acc_req(i) > 0) {
        req_thrust_acc_min(i) = acc_req(i) + G;
        req_thrust_acc_max(i) = max_acc_norm;
      } else {
        req_thrust_acc_max(i) = acc_req(i) + G;
        req_thrust_acc_min(i) = -max_acc_norm;
      }
    } else {
      // x, y
      if (acc_req(i) > 0) {
        req_thrust_acc_min(i) = acc_req(i);
        req_thrust_acc_max(i) = max_acc_norm;
      } else {
        req_thrust_acc_max(i) = acc_req(i);
        req_thrust_acc_min(i) = -max_acc_norm;
      }
    }
  }

  Vector<3> thrust_acc = acc_req - GVEC;
  std::vector<bool> fixed_init{false, false, false};
  reproject_to_sphere(thrust_acc, fixed_init, req_thrust_acc_min,
                      req_thrust_acc_max, acc_req, t_times, max_acc_norm);

  if (fabs(thrust_acc.norm() - max_acc_norm) > 0.01) {
    std::cout << "bad thrust acc norm" << std::endl;
    exit(1);
  }

  Scalar min_tdiff_min_max = std::numeric_limits<Scalar>::max();


  int iter_unimproved_max = 0;
  int last_improved_iter = 0;

  Scalar tmax = std::numeric_limits<Scalar>::max();
  Scalar tmin, tavg;
  Scalar tmax_old = tmax;

  // best trajectory recorded during GD
  std::vector<PMMTrajectory> best_trajectory;
  best_trajectory.resize(3);
  Scalar best_trajectory_tmax = tmax;

  Vector<3> gradients_old = gradients;
  Vector<3> t_times_old = t_times;
  int max_time_idx, max_time_idx_old = 0;
  int min_time_idx, min_time_idx_old = 0;
  Scalar dalph = 0;
  Scalar dalph_old = dalph;
  double default_min_time_change_decay = 1.0;
  double min_time_change_decay = default_min_time_change_decay;
  const double decay_decrease_min_time_swap = 1.0;
  const double decay_decrease_max_time_swap = 0.8;
  const double decay_increase_tmax = 0.9;
  const double decay_decrease_three_constrained = 0.2;

  bool converged = true;
  const int num_iter_opt = NUM_ITERS;
  int iter = 0;
  for (; iter < num_iter_opt; iter++) {
    // std::cout << "------------------" << std::endl;
    // std::cout << "iter " << iter << std::endl;

    converged = true;

    // copy to old vars
    t_times_old = t_times;
    tmax_old = tmax;
    gradients_old = gradients;
    max_time_idx_old = max_time_idx;
    min_time_idx_old = min_time_idx;
    dalph_old = dalph;

    // check if thrust_acc is max
    if (fabs(thrust_acc.norm() - max_acc_norm) > 0.01) {
      std::cout << "bad thrust_dir_size " << thrust_acc.norm() << " vs "
                << max_acc_norm << std::endl;
      exit(1);
    }

    thrust_acc =
      thrust_acc.cwiseMin(req_thrust_acc_max).cwiseMax(req_thrust_acc_min);

    const Vector<3> body_acc = thrust_acc + GVEC;
    const Vector<3> body_acc_down = -thrust_acc + GVEC;
    // acc_req = req_thrust_acc + GVEC;
    for (int i = 0; i < 3; ++i) {
      /*
      if (fabs(body_acc(i)) < fabs(acc_req(i)) &&
          fabs(body_acc(i) - acc_req(i)) < 0.001) {
        body_acc(i) = acc_req(i);
      }
    */

      const Scalar max_acc = body_acc(i);
      const Scalar min_acc = body_acc_down(i);
      const Scalar vs = from.v(i);
      const Scalar ve = to.v(i);
      const Scalar ps = from.p(i);
      const Scalar pe = to.p(i);

      const PMMTrajectory tr(ps, vs, pe, ve, max_acc, min_acc, i, false, true,
                             false);
      // INFO(i << " dtda " << tr.dt_da)

      if (!tr.exists_) {
        std::cout << "iterated to nonexisting one_dim_Scalar_integrator"
                  << std::endl;
        exit(1);
      }

      t_times(i) = tr.time();
      gradients(i) = tr.dt_da_;
      set_axis_trajectory(i, tr);


      // the best a may change sign
      if (i != 2 && acc_req(i) * tr.a_(0) < 0) {
        thrust_acc(i) = tr.a_(0);
        acc_req(i) = copysign(acc_req(i), tr.a_(0));
        const Scalar tmp_acc_max = req_thrust_acc_max(i);
        const Scalar tmp_acc_min = req_thrust_acc_min(i);
        req_thrust_acc_max(i) = -tmp_acc_min;
        req_thrust_acc_min(i) = -tmp_acc_max;
      }
    }

    // use gradient inverse, i.e. da/dt
    Vector<3> gradients_scaled = gradients;  //.cwiseInverse();
    // gradients_scaled.cwiseMax(a_max / 4.0).cwiseMin(-a_max / 4.0);

    tmax = t_times.maxCoeff();

    if (tmax < best_trajectory_tmax) {
      best_trajectory[0] = x_;
      best_trajectory[1] = y_;
      best_trajectory[2] = z_;
      best_trajectory_tmax = tmax;
    }

    get_time_avg_min(tavg, tmin, max_time_idx, min_time_idx, tmax, t_times,
                     gradients_scaled);


    Vector<3> t_times_diff = (t_times_old - t_times);
    t_times_diff = t_times_diff.cwiseAbs().cwiseMax(
      0.0005);  // cap it not tu be too close to 0 for gradient amplification
    Vector<3> gradients_diff =
      (gradients_old - gradients).cwiseAbs().cwiseMax(0.0005);

    bool on_limit = true;

    // setting gradient to zero if on the limit of acc
    while (on_limit) {
      on_limit = false;
      for (int i = 0; i < 3; i++) {
        if (t_times(i) == 0) {
          gradients_scaled(i) = 0;
        }

        if (gradients_scaled(i) != 0) {
          if ((acc_req(i) >= 0 &&
               fabs(thrust_acc(i) - req_thrust_acc_min(i)) < 0.001 &&
               t_times(i) < tavg) ||
              (acc_req(i) <= 0 &&
               fabs(thrust_acc(i) - req_thrust_acc_max(i)) < 0.001 &&
               t_times(i) < tavg)) {
            // on the acc limit while wanting to go further (acc_req different
            // sing from gradient)

            gradients_scaled(i) = 0;

            on_limit = true;
            break;
          }
        }
      }
      get_time_avg_min(tavg, tmin, max_time_idx, min_time_idx, tmax, t_times,
                       gradients_scaled);
    }


    if (iter - last_improved_iter > iter_unimproved_max) {
      iter_unimproved_max = iter - last_improved_iter;
    }
    const Scalar tdiff_min_max = tmax - tmin;
    if (tdiff_min_max < min_tdiff_min_max) {
      min_tdiff_min_max = tdiff_min_max;
      last_improved_iter = iter;
    }

    // end convergence if gradient is zero of >1 axises
    int num_on_limit = 0;
    for (size_t i = 0; i < 3; i++) {
      if (gradients_scaled(i) == 0) {
        num_on_limit++;
      }
    }
    if (num_on_limit > 1) {
      converged = true;
      break;
    }

    // scale the gradient based on time diff to last and signed distance from
    // average
    for (int i = 0; i < 3; i++) {
      Scalar dist_t_avg = fabs(tavg - t_times(i));
      if (t_times(i) > tavg) {
        // INFO("inverse for i " << i)
        gradients_scaled(i) = -gradients_scaled(i) * dist_t_avg /
                              (fabs(t_times_diff(i)) * gradients_diff(i));
      } else {
        // INFO("do not inverse for i " << i)
        gradients_scaled(i) = gradients_scaled(i) * dist_t_avg /
                              (fabs(t_times_diff(i)) * gradients_diff(i));
      }
    }

    // scale alpha angle by distance from average and the dacay
    dalph = min_time_change_decay * (fabs(tmax - tavg) / tmax);
    if (dalph == 0) {
      converged = true;
      break;
    }
    dalph = std::min(0.6, dalph);  // limit the angle

    // so now we have gradients_scaled and need to make them
    // dot(gradients_scaled,thrust_acc_u)==0 and
    // gradient_tang_plane.dot(thrust_acc_u) == 0 &&
    // gradient_tang_plane.norm
    // == thrust acc
    // so we will fix the zeros first

    gradients_scaled = gradients_scaled.normalized() * max_acc_norm;

    // the ones that reached the limit are fixed by force == their scale can
    // not be changed
    std::vector<Scalar> fixed_tang_plane_force = {false, false, false};

    Scalar dot_fixed = 0;
    Scalar g_len_fixed = max_acc_norm * max_acc_norm;
    std::vector<Scalar> fixed_tang_plane = {false, false, false};
    std::vector<int> nf_ids;
    Vector<3> gradient_tang_plane(0, 0, 0);
    for (size_t i = 0; i < 3; i++) {
      if (gradients_scaled(i) == 0) {
        // we have to make acc = acc*cos(alph) + gtp*cos(alph) too keep the
        // acc at the same value
        fixed_tang_plane[i] = true;
        fixed_tang_plane_force[i] = true;
        gradient_tang_plane(i) =
          (thrust_acc(i) - cos(dalph) * thrust_acc(i)) / sin(dalph);

        dot_fixed -= gradient_tang_plane(i) * thrust_acc(i);
        g_len_fixed -= gradient_tang_plane(i) * gradient_tang_plane(i);

      } else {
        nf_ids.push_back(i);
      }
    }

    if (fixed_tang_plane[0] + fixed_tang_plane[1] + fixed_tang_plane[2] == 1) {
      //  now the gradient_tang_plane non fixed parts are set to be
      //  solutions to the
      // system of equations
      // 0 = a_0*g_0 + a_1*g_1 + a_2*g_2
      // a_max^2 = g_0^2 + g_1^2 + g_2^2
      const Scalar acc0_pow2 = (thrust_acc(nf_ids[0]) * thrust_acc(nf_ids[0]));
      const Scalar acc1_pow2 = (thrust_acc(nf_ids[1]) * thrust_acc(nf_ids[1]));
      const Scalar c = -g_len_fixed + (dot_fixed * dot_fixed) / (acc0_pow2);
      const Scalar b = -2.0 * dot_fixed * thrust_acc(nf_ids[1]) / acc0_pow2;
      const Scalar a = (acc1_pow2 + acc0_pow2) / acc0_pow2;
      const Scalar disc = b * b - 4.0 * a * c;
      const Scalar g1_first = (-b + sqrt(disc)) / (2.0 * a);
      const Scalar g1_second = (-b - sqrt(disc)) / (2.0 * a);

      const Scalar g0_first =
        (dot_fixed - thrust_acc(nf_ids[1]) * g1_first) / thrust_acc(nf_ids[0]);
      const Scalar g0_second =
        (dot_fixed - thrust_acc(nf_ids[1]) * g1_second) / thrust_acc(nf_ids[0]);

      // first solution of the quadratic function is the correct one == the
      // one with same sign as gradient scaled


      if (g1_first * gradients_scaled(nf_ids[1]) >= 0 &&
          g0_first * gradients_scaled(nf_ids[0]) >= 0) {
        gradient_tang_plane(nf_ids[1]) = g1_first;
        gradient_tang_plane(nf_ids[0]) = g0_first;
      } else if (g1_second * gradients_scaled(nf_ids[1]) >= 0 &&
                 g0_second * gradients_scaled(nf_ids[0]) >= 0) {
        gradient_tang_plane(nf_ids[1]) = g1_second;
        gradient_tang_plane(nf_ids[0]) = g0_second;
      } else {
        if (max_time_idx == nf_ids[0] &&
            (g0_first * gradients_scaled(nf_ids[0]) >= 0 ||
             g0_second * gradients_scaled(nf_ids[0]) >= 0)) {
          // select the one solution that preserves sign of max_time_idx
          // INFO("here1")
          if (g0_first * gradients_scaled(nf_ids[0]) >= 0) {
            gradient_tang_plane(nf_ids[1]) = g1_first;
            gradient_tang_plane(nf_ids[0]) = g0_first;
          } else if (g0_second * gradients_scaled(nf_ids[0]) >= 0) {
            gradient_tang_plane(nf_ids[1]) = g1_second;
            gradient_tang_plane(nf_ids[0]) = g0_second;
          }
        } else if (max_time_idx == nf_ids[1] &&
                   (g1_first * gradients_scaled(nf_ids[1]) >= 0 ||
                    g1_second * gradients_scaled(nf_ids[1]) >= 0)) {
          // select the one solution that preserves sign of max_time_idx
          // INFO("here2")
          if (g1_first * gradients_scaled(nf_ids[1]) >= 0) {
            gradient_tang_plane(nf_ids[1]) = g1_first;
            gradient_tang_plane(nf_ids[0]) = g0_first;
          } else if (g1_second * gradients_scaled(nf_ids[1]) >= 0) {
            gradient_tang_plane(nf_ids[1]) = g1_second;
            gradient_tang_plane(nf_ids[0]) = g0_second;
          }
        } else if (g1_first * gradients_scaled(nf_ids[1]) >= 0 ||
                   g0_first * gradients_scaled(nf_ids[0]) >= 0) {
          // select any solution that preserves the sign
          // INFO("here3")
          gradient_tang_plane(nf_ids[1]) = g1_first;
          gradient_tang_plane(nf_ids[0]) = g0_first;
        } else if (g1_second * gradients_scaled(nf_ids[1]) >= 0 ||
                   g0_second * gradients_scaled(nf_ids[0]) >= 0) {
          // select any solution that preserves the sign
          // INFO("here4")
          gradient_tang_plane(nf_ids[1]) = g1_second;
          gradient_tang_plane(nf_ids[0]) = g0_second;
        } else {
          std::cout << "there is no max among the nonfixed???" << std::endl;
          std::cout << "non existing same sign gradient plane" << std::endl;

          exit(1);
        }
      }

    } else {
      // fix max first
      gradient_tang_plane =
        gradients_scaled -
        gradients_scaled.normalized().dot(thrust_acc.normalized()) * thrust_acc;
      gradient_tang_plane = gradient_tang_plane.normalized() * max_acc_norm;

      // limit the dalph to be within the scale of gradient to tangent
      // projection
      Scalar angle_between = acos(
        gradient_tang_plane.normalized().dot(gradients_scaled.normalized()));
      if (M_PI_2 - angle_between > 0) {
        dalph = std::min(dalph, M_PI_2 - angle_between);
      }
    }


    Vector<3> new_thrust_acc =
      thrust_acc * cos(dalph) + gradient_tang_plane * sin(dalph);

    // projecting to constrained acc
    std::vector<bool> fixed{false, false, false};
    reproject_to_sphere(new_thrust_acc, fixed, req_thrust_acc_min,
                        req_thrust_acc_max, acc_req, t_times, max_acc_norm);


    const int num_fixed = fixed[0] + fixed[1] + fixed[2];

    if (max_time_idx != max_time_idx_old && iter > 0) {
      // changed max time axis decay
      min_time_change_decay *= decay_decrease_max_time_swap;
    } else if (num_fixed >= 3) {
      // decrease decay num fixed == reached boundary condition
      min_time_change_decay *= decay_decrease_three_constrained;
    } else if (tmax - tmax_old > 0) {
      // decrease decay tmax increased
      min_time_change_decay *= decay_increase_tmax;
    }

    thrust_acc = new_thrust_acc;

    if (fabs(thrust_acc.norm() - max_acc_norm) > 0.01) {
      std::cout << "2!bad thrust_dir_size " << thrust_acc.norm() << " vs "
                << max_acc_norm << std::endl;
      exit(1);
    }


    const Scalar max_time_tst = t_times.maxCoeff();
    const Scalar min_time_tst = t_times.minCoeff();
    if (max_time_tst - min_time_tst > ALLOWED_DIFF_TIMES_RATIO * tmax) {
      converged = false;
    }

    if (converged) {
      break;
    }
  }

  if (tmax > best_trajectory_tmax) {
    // std::cout << "better time during optimization " << tmax << " "
    //           << best_trajectory_tmax << std::endl;
    // exit(1);
    x_ = best_trajectory[0];
    y_ = best_trajectory[1];
    z_ = best_trajectory[2];
  }

  if (equalize_time) {
    const Scalar tr_time = time();
    // std::cout << "equalize time to " << tr_time << std::endl;
    for (size_t i = 0; i < 3; i++) {
      if (get_axis_trajectory(i).time() != tr_time) {
        PMMTrajectory scaled = PMMTrajectory(get_axis_trajectory(i), tr_time);
        set_axis_trajectory(i, scaled);
      }
    }
  } else {
    // std::cout << "not equalizing time " << std::endl;
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

Vector<3> PointMassTrajectory3D::start_acc() const {
  return Vector<3>(x_.a_(0), y_.a_(0), z_.a_(0));
}

Vector<3> PointMassTrajectory3D::end_acc() const {
  return Vector<3>(x_.a_(1), y_.a_(1), z_.a_(1));
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
  const Scalar ax_pow2 = a(0) * a(0);
  const Scalar ay_pow2 = a(1) * a(1);
  const Scalar az_pow2 = a(2) * a(2);
  const Scalar vx_pow2 = v(0) * v(0);
  const Scalar vy_pow2 = v(1) * v(1);
  const Scalar vz_pow2 = v(2) * v(2);

  Scalar t = tto;
  Scalar t_pow2 = t * t;

  Scalar tmp = sqrt(ax_pow2 * t_pow2 + ay_pow2 * t_pow2 + az_pow2 * t_pow2 +
                    2 * a(0) * t * v(0) + vx_pow2 + 2 * a(1) * t * v(1) +
                    vy_pow2 + 2 * a(2) * t * v(2) + vz_pow2);
  Scalar logpart =
    log(ax_pow2 * t + ay_pow2 * t + az_pow2 * t + a(0) * v(0) + a(1) * v(1) +
        a(2) * v(2) + sqrt(ax_pow2 + ay_pow2 + az_pow2) * tmp);
  if (v(0) == 0 && v(1) == 0 && v(2) == 0 && t == 0) {
    logpart = 0;
  }
  const Scalar ds_to =
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

  const Scalar ds_from =
    (0.5) *
    (tmp * (t + (a(0) * v(0) + a(1) * v(1) + a(2) * v(2)) /
                  (ax_pow2 + ay_pow2 + az_pow2)) +
     (1.0 / pow(ax_pow2 + ay_pow2 + az_pow2, 1.5)) *
       ((az_pow2 * (vx_pow2 + vy_pow2) - 2.0 * a(0) * a(2) * v(0) * v(2) -
         2.0 * a(1) * v(1) * (a(0) * v(0) + a(2) * v(2)) +
         ay_pow2 * (vx_pow2 + vz_pow2) + ax_pow2 * (vy_pow2 + vz_pow2)) *
        logpart));

  const Scalar ds = ds_to - ds_from;

  if (!std::isfinite(ds) || ds < 0) {
    std::cout << "non finite ds or bellow zero " << ds << std::endl;
    std::cout << "p " << p.transpose() << std::endl;
    std::cout << "v " << v.transpose() << std::endl;
    std::cout << "a " << a.transpose() << std::endl;
    std::cout << "ds_to " << ds_to << std::endl;
    std::cout << "ds_from " << ds_from << std::endl;
    exit(1);
    return INF;
  }
  return ds;
}

inline void PointMassTrajectory3D::reproject_to_sphere(
  Vector<3> &new_thrust_acc, std::vector<bool> &fixed,
  const Vector<3> &req_thrust_acc_min, const Vector<3> &req_thrust_acc_max,
  const Vector<3> &acc_req, const Vector<3> &t_times, const Scalar &a_max) {
  bool within_constrints = false;
  while (!within_constrints) {
    within_constrints = true;
    if (fabs(new_thrust_acc.norm() - a_max) > 0.0001) {
      within_constrints = false;
    }
    Scalar fixed_len_squared = a_max * a_max;
    Scalar non_fixed_len_squared = 0;
    for (size_t i = 0; i < 3; i++) {
      if (t_times(i) == 0) {
        // t_times(i) == 0 means no acc optimization needed and set to min
        // required one is possible
        if (!fixed[i]) {
          within_constrints = false;
        }
        if (acc_req(i) > 0) {
          new_thrust_acc(i) = req_thrust_acc_min(i);
          fixed[i] = true;
          fixed_len_squared -= req_thrust_acc_min(i) * req_thrust_acc_min(i);
        } else {
          new_thrust_acc(i) = req_thrust_acc_max(i);
          fixed[i] = true;
          fixed_len_squared -= req_thrust_acc_max(i) * req_thrust_acc_max(i);
        }

      } else {
        if (acc_req(i) < 0) {
          if (new_thrust_acc(i) > req_thrust_acc_max(i)) {
            within_constrints = false;
          }
          if (new_thrust_acc(i) > req_thrust_acc_max(i) || fixed[i]) {
            new_thrust_acc(i) = req_thrust_acc_max(i);
            fixed[i] = true;
            fixed_len_squared -= req_thrust_acc_max(i) * req_thrust_acc_max(i);
            // INFO("fixing acc " << i)
          } else {
            non_fixed_len_squared += new_thrust_acc(i) * new_thrust_acc(i);
          }
        } else {
          if (new_thrust_acc(i) < req_thrust_acc_min(i)) {
            within_constrints = false;
          }
          if (new_thrust_acc(i) < req_thrust_acc_min(i) || fixed[i]) {
            new_thrust_acc(i) = req_thrust_acc_min(i);
            fixed[i] = true;
            fixed_len_squared -= req_thrust_acc_min(i) * req_thrust_acc_min(i);
            //  INFO("fixing acc " << i)
          } else {
            non_fixed_len_squared += new_thrust_acc(i) * new_thrust_acc(i);
          }
        }
      }
    }

    if (non_fixed_len_squared != 0) {
      Scalar add_acc = 1.0;
      if (!fixed[2]) {
        non_fixed_len_squared =
          non_fixed_len_squared - 2 * (new_thrust_acc(2) - G) * G - G * G;

        const Scalar c = -fixed_len_squared + G * G;
        const Scalar b = 2 * (new_thrust_acc(2) - G) * G;
        const Scalar a = non_fixed_len_squared;
        add_acc = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
      } else {
        add_acc = sqrt(fixed_len_squared / non_fixed_len_squared);
      }

      for (size_t i = 0; i < 3; i++) {
        if (!fixed[i]) {
          if (i == 2) {
            new_thrust_acc(i) = (new_thrust_acc(i) - G) * add_acc + G;
          } else {
            new_thrust_acc(i) = new_thrust_acc(i) * add_acc;
          }
        }
      }
    } else {
      Scalar scale_non_fixed = 1.0;
      scale_non_fixed =
        sqrt((a_max * a_max) / ((a_max * a_max) - fixed_len_squared));
      // INFO("2scale_non_fixed " << scale_non_fixed)
      for (size_t i = 0; i < 3; i++) {
        new_thrust_acc(i) = new_thrust_acc(i) * (scale_non_fixed);
      }
    }
  }
}

inline void PointMassTrajectory3D::get_time_avg_min(
  Scalar &tavg, Scalar &tmin, int &max_time_idx, int &min_time_idx,
  const Scalar &tmax, const Vector<3> &t_times,
  const Vector<3> &gradients_scaled) {
  tmin = DBL_MAX;
  tavg = 0;
  Scalar tavg_num = 0;
  for (size_t i = 0; i < 3; i++) {
    if (t_times(i) != 0 && gradients_scaled(i) != 0) {
      tavg += t_times(i);
      tavg_num += 1.0;
      if (t_times(i) < tmin) {
        tmin = t_times(i);
        min_time_idx = i;
      }
      if (t_times(i) == tmax) {
        max_time_idx = i;
      }
    }
  }
  tavg /= tavg_num;
  if (max_time_idx == -1 || min_time_idx == -1) {
    std::cout << "wrong max_time_idx or min_time_idx " << max_time_idx << " "
              << min_time_idx << std::endl;
    exit(1);
  }
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
      std::cout << "bad axis index " << i << std::endl;
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
      std::cout << "bad axis index " << i << std::endl;
      return x_;
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