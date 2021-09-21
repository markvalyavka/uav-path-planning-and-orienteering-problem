#include "pmm_trajectory3d.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace agi {

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             const Vector<3> max_acc,
                                             const bool equalize_time) {
  x_ = PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc(0),
                     -max_acc(0), 0);
  y_ = PMMTrajectory(from.p(1), from.v(1), to.p(1), to.v(1), max_acc(1),
                     -max_acc(1), 1);
  z_ = PMMTrajectory(from.p(2), from.v(2), to.p(2), to.v(2), max_acc(2),
                     -max_acc(2), 2);

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
                                             Scalar max_acc_norm,
                                             const bool equalize_time) {
  static const Scalar ALLOWED_DIFF_TIMES_RATIO{0.0001};
  static const Scalar NUM_ITERS{30};

  Vector<3> t_times(0.0, 0.0, 0.0);
  Vector<3> gradients(0.0, 0.0, 0.0);
  Vector<3> acc_req(0.01, 0.01, 0.01);


  double tmax = std::numeric_limits<double>::max();

  Vector<3> t_times_old = t_times;
  Vector<3> gradients_old = gradients;
  int max_time_idx, max_time_idx_old = 0;
  int min_time_idx, min_time_idx_old = 0;
  double tmax_old = tmax;
  double dalph = 0;
  double dalph_old = dalph;

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
    if (fabs(thrust_acc.norm() - a_max) > 0.01) {
      std::cout << "bad thrust_dir_size " << thrust_acc.norm() << " vs "
                << a_max << std::endl;
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

      const double max_acc = body_acc(i);
      const double min_acc = body_acc_down(i);
      const double vs = from_state.v(i);
      const double ve = to_state.v(i);
      const double ps = from_state.p(i);
      const double pe = to_state.p(i);

      const Tr1D tr =
        one_dim_double_integrator_two_acc(ps, vs, pe, ve, max_acc, min_acc, i);
      // INFO(i << " dtda " << tr.dt_da)

      if (!tr.exists) {
        ERROR_RED("iterated to nonexisting one_dim_double_integrator");
        INFO("axis " << i)
        INFO_VAR(tr.a1)
        INFO_VAR(tr.a2)
        INFO_VAR(max_acc)

        INFO_VAR(req_thrust_acc_max.transpose())
        INFO_VAR(req_thrust_acc_min.transpose())
        INFO_VAR(acc_req.transpose())
        INFO_VAR(acc_req.norm())
        INFO_VAR(body_acc.transpose())
        INFO_VAR(thrust_acc.transpose())
        INFO_VAR(thrust_acc.norm())
        INFO_VAR(a_max)
        INFO_VAR(acc_req.transpose())
        exit(1);
      }

      t_times(i) = tr.time();
      gradients(i) = tr.dt_da;
      tr_max_acc.set_by_axis(i, tr);

      // if (i == 2) {
      //   gradients(i) = copysign(gradients(i), -1);
      // }

      // the best a may change sign
      if (i != 2 && acc_req(i) * tr.a1 < 0) {
        thrust_acc(i) = tr.a1;
        acc_req(i) = copysign(acc_req(i), tr.a1);
        const double tmp_acc_max = req_thrust_acc_max(i);
        const double tmp_acc_min = req_thrust_acc_min(i);
        req_thrust_acc_max(i) = -tmp_acc_min;
        req_thrust_acc_min(i) = -tmp_acc_max;
      }
    }

    // use gradient inverse, i.e. da/dt
    Vector<3> gradients_scaled = gradients;  //.cwiseInverse();
    // gradients_scaled.cwiseMax(a_max / 4.0).cwiseMin(-a_max / 4.0);

    tmax = t_times.maxCoeff();

    get_time_avg_min(tavg, tmin, max_time_idx, min_time_idx, tmax, t_times,
                     gradients_scaled);


    INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_CYAN, "tmax " << tmax)
    INFO_COND_COLOR(
      PRINT_DEBUG, OUTPUT_CYAN,
      "t_times " << t_times(0) << " " << t_times(1) << " " << t_times(2))
    INFO_COND(PRINT_DEBUG, "gradients " << gradients(0) << " " << gradients(1)
                                        << " " << gradients(2))
    INFO_COND(PRINT_DEBUG, "gradients_scaled init "
                             << gradients_scaled(0) << " "
                             << gradients_scaled(1) << " "
                             << gradients_scaled(2))
    INFO_COND(PRINT_DEBUG, "thrust_acc " << thrust_acc(0) << " "
                                         << thrust_acc(1) << " "
                                         << thrust_acc(2))
    INFO_COND(PRINT_DEBUG, "tmin " << tmin)


    Vector<3> t_times_diff = (t_times_old - t_times);
    t_times_diff = t_times_diff.cwiseAbs().cwiseMax(
      0.0005);  // cap it not tu be too close to 0 for gradient amplification
    Vector<3> gradients_diff =
      (gradients_old - gradients).cwiseAbs().cwiseMax(0.0005);

    INFO_COND(PRINT_DEBUG, "t_times_diff " << t_times_diff.transpose())
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

            INFO_COND_COLOR(
              PRINT_DEBUG, OUTPUT_RED,
              "is on limit of descent, do not consider for  avg.... " << i)
            INFO_COND(PRINT_DEBUG, "acc_req " << acc_req.transpose())


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
    const double tdiff_min_max = tmax - tmin;
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
      double dist_t_avg = fabs(tavg - t_times(i));
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


    INFO_COND(PRINT_DEBUG, "gradients_scaled bef " << gradients_scaled(0) << " "
                                                   << gradients_scaled(1) << " "
                                                   << gradients_scaled(2))
    INFO_COND(PRINT_DEBUG,
              "tmax " << tmax << " tavg " << tavg << " tmin " << tmin)

    // scale alpha angle by distance from average and the dacay
    dalph = min_time_change_decay * (fabs(tmax - tavg) / tmax);
    if (dalph == 0) {
      converged = true;
      break;
    }
    dalph = std::min(0.6, dalph);  // limit the angle
    INFO_COND(PRINT_DEBUG, "dalph " << dalph)

    // so now we have gradients_scaled and need to make them
    // dot(gradients_scaled,thrust_acc_u)==0 and
    // gradient_tang_plane.dot(thrust_acc_u) == 0 &&
    // gradient_tang_plane.norm
    // == thrust acc
    // so we will fix the zeros first

    gradients_scaled = gradients_scaled.normalized() * a_max;

    INFO_COND(PRINT_DEBUG, "gradients_scaled normalized "
                             << gradients_scaled(0) << " "
                             << gradients_scaled(1) << " "
                             << gradients_scaled(2))

    // the ones that reached the limit are fixed by force == their scale can
    // not be changed
    std::vector<double> fixed_tang_plane_force = {false, false, false};

    double dot_fixed = 0;
    double g_len_fixed = a_max * a_max;
    std::vector<double> fixed_tang_plane = {false, false, false};
    std::vector<int> nf_ids;
    Vector<3> gradient_tang_plane(0, 0, 0);
    for (size_t i = 0; i < 3; i++) {
      if (gradients_scaled(i) == 0) {
        // we have to make acc = acc*cos(alph) + gtp*cos(alph) too keep the
        // acc at the same value
        INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_MAGENTA, "on limit for idx " << i)
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
      const double acc0_pow2 = (thrust_acc(nf_ids[0]) * thrust_acc(nf_ids[0]));
      const double acc1_pow2 = (thrust_acc(nf_ids[1]) * thrust_acc(nf_ids[1]));
      const double c = -g_len_fixed + (dot_fixed * dot_fixed) / (acc0_pow2);
      const double b = -2.0 * dot_fixed * thrust_acc(nf_ids[1]) / acc0_pow2;
      const double a = (acc1_pow2 + acc0_pow2) / acc0_pow2;
      const double disc = b * b - 4.0 * a * c;
      const double g1_first = (-b + sqrt(disc)) / (2.0 * a);
      const double g1_second = (-b - sqrt(disc)) / (2.0 * a);

      const double g0_first =
        (dot_fixed - thrust_acc(nf_ids[1]) * g1_first) / thrust_acc(nf_ids[0]);
      const double g0_second =
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
          INFO("there is no max among the nonfixed???")
          INFO("gradients_scaled normalized " << gradients_scaled(0) << " "
                                              << gradients_scaled(1) << " "
                                              << gradients_scaled(2))
          INFO("gradients " << gradients(0) << " " << gradients(1) << " "
                            << gradients(2))
          INFO("non existing same sign gradient plane")
          INFO_VAR(g1_first)
          INFO_VAR(g1_second)
          INFO_VAR(gradients_scaled(nf_ids[1]))

          INFO_VAR(tavg);
          INFO_VAR(g0_first)
          INFO_VAR(g0_second)
          INFO_VAR(gradients_scaled(nf_ids[0]))
          INFO_VAR(tavg);
          INFO_VAR(disc)

          INFO_VAR(req_thrust_acc_max.transpose())
          INFO_VAR(req_thrust_acc_min.transpose())

          INFO_VAR(req_max_thrust_acc.transpose())
          INFO_VAR(req_max_thrust_acc.norm())
          INFO_VAR(a_max)

          INFO("t_times " << t_times(0) << " " << t_times(1) << " "
                          << t_times(2))
          INFO("t_times-tavg " << t_times(0) - tavg << " " << t_times(1) - tavg
                               << " " << t_times(2) - tavg)
          INFO("t_times_diff " << t_times_diff.transpose())
          INFO("thrust_acc " << thrust_acc(0) << " " << thrust_acc(1) << " "
                             << thrust_acc(2) << " norm " << thrust_acc.norm()
                             << " a_max " << a_max)
          INFO_VAR(nf_ids[0])
          INFO_VAR(nf_ids[1])
          INFO_VAR(max_time_idx)
          INFO_VAR(min_time_idx)
          exit(1);
        }
      }


    } else {
      // fix max first
      // INFO("thrust_acc norm " << thrust_acc.norm())
      gradient_tang_plane =
        gradients_scaled -
        gradients_scaled.normalized().dot(thrust_acc.normalized()) * thrust_acc;
      gradient_tang_plane = gradient_tang_plane.normalized() * a_max;

      // limit the dalph to be within the scale of gradient to tangent
      // projection
      double angle_between = acos(
        gradient_tang_plane.normalized().dot(gradients_scaled.normalized()));
      if (M_PI_2 - angle_between > 0) {
        dalph = std::min(dalph, M_PI_2 - angle_between);
        INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_RED,
                        "!!!!!!!! limiting dalph to " << dalph)
      }
      // INFO_MAGENTA("angle_between " << angle_between)
    }

    INFO_COND(PRINT_DEBUG, "gradient_tang_plane dot_equal "
                             << gradient_tang_plane.transpose() << " norm "
                             << gradient_tang_plane.norm())


    INFO_COND(PRINT_DEBUG, "acc_req " << acc_req(0) << " " << acc_req(1) << " "
                                      << acc_req(2))


    INFO_COND(PRINT_DEBUG, "gradient_tang_plane "
                             << gradient_tang_plane.transpose() << " norm "
                             << gradient_tang_plane.norm())

    Vector<3> new_thrust_acc =
      thrust_acc * cos(dalph) + gradient_tang_plane * sin(dalph);


    INFO_COND(PRINT_DEBUG, "new_thrust_acc new "
                             << new_thrust_acc(0) << " " << new_thrust_acc(1)
                             << " " << new_thrust_acc(2) << " norm "
                             << new_thrust_acc.norm() << " a_max " << a_max)


    // INFO("thrust_acc_tst size " << new_thrust_acc.norm() << " a_max "
    //                             << a_max)

    // projecting to constrained acc
    std::vector<bool> fixed{false, false, false};
    reproject_to_sphere(new_thrust_acc, fixed, req_thrust_acc_min,
                        req_thrust_acc_max, acc_req, t_times, a_max);


    const int num_fixed = fixed[0] + fixed[1] + fixed[2];

    if (max_time_idx != max_time_idx_old && iter > 0) {
      min_time_change_decay *= decay_decrease_max_time_swap;
      INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_GREEN,
                      "changed max time axis decay " << min_time_change_decay)
    } else if (num_fixed >= 3) {
      min_time_change_decay *= decay_decrease_three_constrained;
      INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_GREEN,
                      "decrease decay num fixed " << min_time_change_decay)
    } else if (tmax - tmax_old > 0) {
      min_time_change_decay *= decay_increase_tmax;
      INFO_COND_COLOR(PRINT_DEBUG, OUTPUT_GREEN,
                      "decrease decay tmax increased " << min_time_change_decay)
    }


    thrust_acc = new_thrust_acc;

    if (fabs(thrust_acc.norm() - a_max) > 0.01) {
      INFO("2!bad thrust_dir_size " << thrust_acc.norm() << " vs " << a_max);
      INFO_VAR(iter)
      INFO_VAR(thrust_acc.norm())
      INFO_VAR(a_max)
      exit(1);
    }


    const double max_time_tst = t_times.maxCoeff();
    const double min_time_tst = t_times.minCoeff();
    if (max_time_tst - min_time_tst > ALLOWED_DIFF_TIMES_RATIO * tmax) {
      converged = false;
    }

    if (converged) {
      break;
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
  // const Scalar &sx = p(0);
  // const Scalar &sy = p(1);
  // const Scalar &sz = p(2);

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