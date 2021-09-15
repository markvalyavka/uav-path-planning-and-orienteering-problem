
#include "velocity_search_graph.hpp"

#include <cfloat>

#include "agilib/utils/timer.hpp"

namespace agi {

VelocitySearchGraph::VelocitySearchGraph(
  const Scalar max_acc, const Scalar max_yaw_pitch_ang,
  const Scalar precision_yaw_pitch_ang,
  const Scalar yaw_pitch_cone_angle_boundary,
  const Scalar min_velocity_size_boundary, const Scalar min_velocity_size,
  const Scalar max_velocity_size, const Scalar precision_velocity_size)
  : max_acc_(max_acc),
    max_yaw_pitch_ang_(max_yaw_pitch_ang),
    precision_yaw_pitch_ang_(precision_yaw_pitch_ang),
    yaw_pitch_cone_angle_boundary_(yaw_pitch_cone_angle_boundary),
    min_velocity_size_boundary_(min_velocity_size_boundary),
    min_velocity_size_(min_velocity_size),
    max_velocity_size_(max_velocity_size),
    precision_velocity_size_(precision_velocity_size) {}


MultiWaypointTrajectory VelocitySearchGraph::find_velocities_in_positions(
  const std::vector<Vector<3>>& gates_waypoints,
  const Vector<3>& start_velocity, const Vector<3>& end_velocity,
  const std::vector<Scalar>& gates_yaw_deg, const bool end_free) {
  const int gates_size = gates_waypoints.size();

  std::vector<Vector<3>> found_gates_speeds;
  std::vector<Scalar> found_gates_times;
  // also allow optimizing the end but omit the start
  found_gates_speeds.resize(gates_size - 1);
  found_gates_times.resize(gates_size - 1);

  std::vector<std::vector<Vector<2>>> gates_yaw_pitch_size_ranges;
  gates_yaw_pitch_size_ranges.resize(
    gates_size - 1, {Vector<2>(-max_yaw_pitch_ang_, max_yaw_pitch_ang_),
                     Vector<2>(-max_yaw_pitch_ang_, max_yaw_pitch_ang_),
                     Vector<2>(min_velocity_size_, max_velocity_size_)});


  // samles have [gateid][sampleid][{velocity,(yaw, pitch,
  // vel_size),(yaw_coneindex, pitch_coneindex, vel_size_coneindex)}]
  std::vector<std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>>
    gate_velocity_samples;
  gate_velocity_samples.resize(gates_size);
  gate_velocity_samples[0].push_back(
    {start_velocity, Vector<3>::Zero(), Vector<3>::Zero()});

  // shortest time structure for samples, add also the start, pair ==
  // {id_of_sample_in_gate_before,time from the start}
  std::vector<std::vector<std::pair<int, Scalar>>> shortest_samples_times;
  shortest_samples_times.resize(gate_velocity_samples.size());
  shortest_samples_times[0].push_back({-1, 0});

  const int num_vel_samples_yaw_pitch_ang =
    ceil(2 * max_yaw_pitch_ang_ / precision_yaw_pitch_ang_) +
    1;  // 2* range+zero
  const Scalar vel_samples_yaw_pitch_ang_step =
    num_vel_samples_yaw_pitch_ang > 1
      ? 2 * max_yaw_pitch_ang_ / ((Scalar)(num_vel_samples_yaw_pitch_ang - 1))
      : 0;

  const int num_vel_samples_size =
    ceil((max_velocity_size_ - min_velocity_size_) / precision_velocity_size_) +
    1;  // from 0 to max velocity size

  const Scalar vel_samples_size_step =
    num_vel_samples_size > 1 ? (max_velocity_size_ - min_velocity_size_) /
                                 ((Scalar)(num_vel_samples_size - 1))
                             : 0;

  int sample_num_gates = gates_size - 2;
  if (end_free) {
    sample_num_gates += 1;
  } else {
    gate_velocity_samples[gates_size - 1].push_back(
      {end_velocity, Vector<3>::Zero(), Vector<3>::Zero()});
    shortest_samples_times[gate_velocity_samples.size() - 1].push_back(
      {-1, DBL_MAX});  // add end
  }


  Timer timer;
  timer.tic();


  const int num_iter_improve = 5;

  // iteration over improvements of velocities == focusing velocity cones
  for (size_t iterimp = 0; iterimp < num_iter_improve; iterimp++) {
    // loop gates 1....sample_num_gates (0 is start and gatesize-1 is end)
    // here creating only the samples
    for (size_t gid = 1; gid <= sample_num_gates; gid++) {
      const int num_samples_gate =
        num_vel_samples_yaw_pitch_ang * num_vel_samples_yaw_pitch_ang *
        num_vel_samples_size;  // zero size does not need yaw-pitch samples

      // get current min-max angles(relative to gate yaw) and velocity sizes
      const Vector<2> min_max_yaw = gates_yaw_pitch_size_ranges[gid - 1][0];
      const Vector<2> min_max_pitch = gates_yaw_pitch_size_ranges[gid - 1][1];
      const Vector<2> min_max_size = gates_yaw_pitch_size_ranges[gid - 1][2];

      // angle step sizes
      const Scalar gate_vel_samples_size_step =
        (min_max_size(1) - min_max_size(0)) /
        ((Scalar)(num_vel_samples_size - 1));
      const Scalar gate_vel_samples_yaw_step =
        (min_max_yaw(1) - min_max_yaw(0)) /
        ((Scalar)(num_vel_samples_yaw_pitch_ang - 1));
      const Scalar gate_vel_samples_pitch_step =
        (min_max_pitch(1) - min_max_pitch(0)) /
        ((Scalar)(num_vel_samples_yaw_pitch_ang - 1));

      const Scalar gate_yaw = gates_yaw_deg[gid];

      gate_velocity_samples[gid].resize(num_samples_gate);

      shortest_samples_times[gid].resize(num_samples_gate, {-1, DBL_MAX});
      const std::pair<int, Scalar> not_reached = std::make_pair(-1, DBL_MAX);
      std::fill(shortest_samples_times[gid].begin(),
                shortest_samples_times[gid].end(), not_reached);

      // save samples values to the gate_velocity_samples
      int i = 0;
      for (size_t z_s = 0; z_s < num_vel_samples_size; z_s++) {
        const Scalar vel_size =
          min_max_size(0) + z_s * gate_vel_samples_size_step;
        for (size_t x_s = 0; x_s < num_vel_samples_yaw_pitch_ang; x_s++) {
          const Scalar yaw = min_max_yaw(0) + x_s * gate_vel_samples_yaw_step;
          for (size_t y_s = 0; y_s < num_vel_samples_yaw_pitch_ang; y_s++) {
            const Scalar pitch =
              min_max_pitch(0) + y_s * gate_vel_samples_pitch_step;
            Matrix<3, 3> m;
            m = Eigen::AngleAxisd(pitch * M_PI / 180.0, Vector<3>::UnitY()) *
                Eigen::AngleAxisd((gate_yaw + yaw) * M_PI / 180.0,
                                  Vector<3>::UnitZ());
            const Vector<3> vec = m * Vector<3>::UnitX() * vel_size;
            const Vector<3> vec_yaw_pitch_size(yaw, pitch, vel_size);
            const Vector<3> vec_yaw_pitch_size_idx(x_s, y_s, z_s);

            std::get<0>(gate_velocity_samples[gid][i]) = vec;
            std::get<1>(gate_velocity_samples[gid][i]) = vec_yaw_pitch_size;
            std::get<2>(gate_velocity_samples[gid][i]) = vec_yaw_pitch_size_idx;
            i++;
          }
        }
      }
    }

    // now we have the samples, search the shortest time path through the graph
    // loop from the first gate gid_to = 1
    for (size_t gid_to = 1; gid_to < gates_size; gid_to++) {
      // get locally the parts
      const std::vector<std::pair<int, Scalar>>& from_shortest_samples_dists =
        shortest_samples_times[gid_to - 1];
      std::vector<std::pair<int, Scalar>>& to_shortest_samples_dists =
        shortest_samples_times[gid_to];
      const std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>&
        from_samples = gate_velocity_samples[gid_to - 1];
      const std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>&
        to_samples = gate_velocity_samples[gid_to];

      // get the positions of the gates/start/end
      const Vector<3>& from_p = gates_waypoints[gid_to - 1];
      const Vector<3>& to_p = gates_waypoints[gid_to];

      Scalar min_dist_to_gate = DBL_MAX;

      // loop for indexes of the samples in the gate (gid_to-1)
      for (int idx_from = 0; idx_from < from_shortest_samples_dists.size();
           idx_from++) {
        const Scalar& time_from = from_shortest_samples_dists[idx_from].second;
        const Vector<3>& from_v = std::get<0>(from_samples[idx_from]);
        if (time_from == DBL_MAX) {
          // no path reaching that sample
          continue;
        }

        // loop for indexes of the samples in the gate gid_to
        for (int idx_to = 0; idx_to < to_shortest_samples_dists.size();
             idx_to++) {
          std::pair<int, Scalar>& shortest_to =
            to_shortest_samples_dists[idx_to];

          // get the velocitions in particular samples
          const Vector<3>& to_v = std::get<0>(to_samples[idx_to]);

          QuadState from_state;
          from_state.p = from_p;
          from_state.v = from_v;
          QuadState to_state;
          to_state.p = to_p;
          to_state.v = to_v;

          const PointMassTrajectory3D tr_max_acc(from_state, to_state, max_acc_,
                                                 false);

          if (tr_max_acc.exists()) {
            const Scalar time_between = tr_max_acc.time();
            const Scalar time_tot = time_between + time_from;

            if (time_tot < min_dist_to_gate) {
              min_dist_to_gate = time_tot;
            }
            if (time_tot < shortest_to.second) {
              shortest_to.first = idx_from;
              shortest_to.second = time_tot;
            }
          }
        }
      }


      if (min_dist_to_gate == DBL_MAX) {
        std::cout << "there is no connection to gate " << gid_to
                  << " with position " << gates_waypoints[gid_to].transpose()
                  << std::endl;
        std::cout << from_p.transpose() << std::endl;
        std::cout << to_p.transpose() << std::endl;
        return MultiWaypointTrajectory();
      }
    }

    // find the shortest-time trajectory among the end samples
    Scalar shortest_time = DBL_MAX;
    int end_best_idx = -1;
    const int endi = shortest_samples_times.size() - 1;
    for (size_t i = 0; i < shortest_samples_times[endi].size(); i++) {
      const Scalar time = shortest_samples_times[endi][i].second;
      if (time < shortest_time) {
        shortest_time = time;
        end_best_idx = i;
      }
    }

    // trace back the trajectory from the shortest final sample
    // additionally also optimize the ranges of the velocity cone
    int prev_sample_idx = end_best_idx;
    for (int g_id = endi; g_id > 0; g_id--) {
      found_gates_speeds[g_id - 1] =
        std::get<0>(gate_velocity_samples[g_id][prev_sample_idx]);

      Vector<3> indexes =
        std::get<2>(gate_velocity_samples[g_id][prev_sample_idx]);

      for (size_t i = 0; i < 3; i++) {
        const Vector<2> current_range =
          gates_yaw_pitch_size_ranges[g_id - 1][i];

        const Scalar range_size_half =
          (current_range(1) - current_range(0)) / 2.0;
        if (indexes(i) == 0) {
          // move the range of minmax down
          gates_yaw_pitch_size_ranges[g_id - 1][i] =
            Vector<2>(current_range(0) - range_size_half,
                      current_range(1) - range_size_half);
          if (i == 2) {
            // limit velocity size to be bigger than min_velocity_size_boundary_
            gates_yaw_pitch_size_ranges[g_id - 1][i] =
              gates_yaw_pitch_size_ranges[g_id - 1][i].cwiseMax(
                min_velocity_size_boundary_);
          } else {
            // limit yaw pitch change from gate direction to be within
            // -yaw_pitch_cone_angle_boundary_ and
            // yaw_pitch_cone_angle_boundary_
            gates_yaw_pitch_size_ranges[g_id - 1][i] =
              gates_yaw_pitch_size_ranges[g_id - 1][i]
                .cwiseMax(-yaw_pitch_cone_angle_boundary_)
                .cwiseMin(yaw_pitch_cone_angle_boundary_);
          }
        } else if (indexes(i) == 2) {
          // move the range of minmax up
          gates_yaw_pitch_size_ranges[g_id - 1][i] =
            Vector<2>(current_range(0) + range_size_half,
                      current_range(1) + range_size_half);
          if (i < 2) {
            // limit yaw pitch change from gate direction to be within
            // -yaw_pitch_cone_angle_boundary_ and
            // yaw_pitch_cone_angle_boundary_
            gates_yaw_pitch_size_ranges[g_id - 1][i] =
              gates_yaw_pitch_size_ranges[g_id - 1][i]
                .cwiseMax(-yaw_pitch_cone_angle_boundary_)
                .cwiseMin(yaw_pitch_cone_angle_boundary_);
          }
        } else {
          // make smaller range around current sample
          gates_yaw_pitch_size_ranges[g_id - 1][i] =
            Vector<2>(current_range(0) + range_size_half / 2.0,
                      current_range(1) - range_size_half / 2.0);
        }
      }

      // traceback previous sample index here
      prev_sample_idx = shortest_samples_times[g_id][prev_sample_idx].first;
    }
  }

  timer.toc();

  // find the shortest among the end ones
  Scalar shortest_time = DBL_MAX;
  int end_best_idx = -1;  // end idx is 0
  const int endi = shortest_samples_times.size() - 1;
  for (size_t i = 0; i < shortest_samples_times[endi].size(); i++) {
    const Scalar time = shortest_samples_times[endi][i].second;
    if (time < shortest_time) {
      shortest_time = time;
      end_best_idx = i;
    }
  }

  int prev_sample_idx = end_best_idx;
  for (int g_id = endi; g_id > 0; g_id--) {
    found_gates_speeds[g_id - 1] =
      std::get<0>(gate_velocity_samples[g_id][prev_sample_idx]);
    found_gates_times[g_id - 1] =
      shortest_samples_times[g_id][prev_sample_idx].second;

    prev_sample_idx = shortest_samples_times[g_id][prev_sample_idx].first;
  }

  MultiWaypointTrajectory trajectories;

  Scalar time_sum = 0;

  QuadState from_start_state;
  from_start_state.p = gates_waypoints[0];
  from_start_state.v = start_velocity;
  QuadState to_state;
  to_state.p = gates_waypoints[1];
  to_state.v = found_gates_speeds[0];
  PointMassTrajectory3D tr_max_acc_start =
    PointMassTrajectory3D(from_start_state, to_state, max_acc_);
  trajectories.push_back(tr_max_acc_start);
  time_sum += tr_max_acc_start.time();

  if (!tr_max_acc_start.exists()) {
    std::cout << "pmm trajectory does not exist from start" << std::endl;
    return MultiWaypointTrajectory();
  }

  for (size_t i = 2; i < gates_waypoints.size(); i++) {
    QuadState from_state_b;
    from_state_b.p = gates_waypoints[i - 1];
    from_state_b.v = found_gates_speeds[i - 2];
    QuadState to_state_b;
    to_state_b.p = gates_waypoints[i];
    to_state_b.v = found_gates_speeds[i - 1];
    PointMassTrajectory3D tr_max_between =
      PointMassTrajectory3D(from_state_b, to_state_b, max_acc_);
    trajectories.push_back(tr_max_between);

    time_sum += tr_max_between.time();
    if (!tr_max_between.exists()) {
      // this trajectory does not exists
      std::cout << "found pmm trajectory does not exist" << std::endl;
      return MultiWaypointTrajectory();
    }
  }


  if (fabs(time_sum - shortest_time) > 0.0001) {
    std::cout << "error calculationg pmm trajectory" << std::endl;
    return MultiWaypointTrajectory();
  }

  return trajectories;
}

void VelocitySearchGraph::save_track_trajectory_equidistant(
  const MultiWaypointTrajectory& trajectories, std::string filename) {
  std::vector<std::vector<QuadState>> samples;
  samples.resize(1);
  const Scalar desired_ds = 0.250;
  Scalar trajectory_length = 0;

  // get the trajectories length
  for (size_t i = 0; i < trajectories.size(); i++) {
    const Scalar ttime = trajectories[i].time_min();
    const Scalar traj_length =
      trajectories[i].get_length_between_times(0, ttime);
    trajectory_length += traj_length;
  }

  int num_samples = ceil(trajectory_length / desired_ds);
  Scalar ds = trajectory_length / ((Scalar)(num_samples - 1));

  QuadState beginstate = trajectories[0].state_in_time(0);
  beginstate.t = 0;
  samples[0].push_back(beginstate);

  int tr_id = 0;
  Scalar tr_start_time = 0;
  Scalar tr_end_time = trajectories[0].time_min();

  Scalar lower_time = 0;
  Scalar upper_time = 0;
  Scalar dist_first = 0;
  Scalar time_middle = 0;
  Scalar t_current = 0;
  Scalar t_from_start = 0;
  Scalar accumulated_s = 0;
  for (size_t i = 1; i <= num_samples; i++) {
    const Scalar s = i * ds;
    Scalar target_ds = s - accumulated_s;

    Scalar dist_to_end =
      trajectories[tr_id].get_length_between_times(t_current, tr_end_time);

    if (dist_to_end < target_ds) {
      tr_id++;
      if (tr_id >= trajectories.size()) {
        break;
      }
      t_current = 0;
      t_from_start += tr_end_time;
      tr_end_time = trajectories[tr_id].time_min();
      target_ds -= dist_to_end;
      accumulated_s += dist_to_end;
    }

    lower_time = t_current;
    upper_time = tr_end_time;

    do {
      time_middle = (lower_time + upper_time) / 2.0;
      dist_first =
        trajectories[tr_id].get_length_between_times(t_current, time_middle);

      if (dist_first > target_ds) {
        // soulution between lower_time and time_middle
        upper_time = time_middle;
      } else {
        // soulution between time_middle and upper_time
        lower_time = time_middle;
      }
    } while (fabs(dist_first - target_ds) > 0.0000000001);

    // add also samples of acc switch
    std::vector<Scalar> acc_switch_times;
    for (size_t axi = 0; axi < 3; axi++) {
      const Scalar t1 = trajectories[tr_id].get_axis_switch_time(axi);
      if (t1 > t_current && t1 < time_middle) {
        acc_switch_times.push_back(t1);
      }
    }
    std::sort(acc_switch_times.begin(), acc_switch_times.end());

    // converged
    t_current = time_middle;
    accumulated_s += dist_first;

    if (t_current < 0 || t_current > trajectories[tr_id].time_min()) {
      std::cout << "bad time" << std::endl;
      return;
    }

    QuadState dronestate = trajectories[tr_id].state_in_time(t_current);
    if (fabs(dronestate.p(1)) > 30) {
      std::cout << "bad distance" << std::endl;
      return;
    }
    dronestate.t = t_from_start + t_current;

    const Scalar dist_to_last = (samples[0].back().p - dronestate.p).norm();
    if (fabs(dist_to_last - ds) > 0.01) {
      std::cout << "bad distance change in save_track_trajectory_equidistant "
                << dist_to_last << std::endl;
      std::cout << "sample " << i << " out of " << num_samples << std::endl;
      return;
    }
    samples[0].push_back(dronestate);
  }

  saveSamplesToFile(filename, samples);
}

void VelocitySearchGraph::save_track_trajectory(
  const MultiWaypointTrajectory& trajectories, std::string filename) {
  std::vector<std::vector<QuadState>> samples;
  samples.resize(1);
  const Scalar dt_desired = 0.01;
  Scalar trajectories_time = 0;
  for (size_t i = 0; i < trajectories.size(); i++)
    trajectories_time += trajectories[i].time();

  const int num_samples = ceil(trajectories_time / dt_desired);

  samples[0].resize(num_samples + 1);
  const Scalar dt = trajectories_time / ((Scalar)(num_samples));
  int tr_id = 0;
  Scalar tr_start_time = 0;
  Scalar tr_end_time = trajectories[0].time();
  for (size_t i = 0; i <= num_samples; i++) {
    const Scalar t = i * dt;
    if (t > tr_end_time) {
      tr_id++;
      if (tr_id >= trajectories.size()) {
        break;
      }
      tr_start_time = tr_end_time;
      tr_end_time = tr_end_time + trajectories[tr_id].time();
    }

    QuadState ds = trajectories[tr_id].state_in_time(t - tr_start_time);
    ds.t = t;
    samples[0][i] = ds;
  }

  saveSamplesToFile(filename, samples);
}

std::vector<Vector<3>> VelSearchGraph::sampleTrajectory(
  const TrMaxAcc3D& trajectory, const Scalar ds_desired) {
  std::vector<Vector<3>> samples;
  // INFO("sample trajectory begin")


  // get the trajectories time and length
  const Scalar ttime = trajectory.timemin();
  // INFO_VAR(ttime)
  const Scalar trajectory_length = trajectory.get_length(0, ttime);
  // INFO_VAR(trajectory_length)

  // define num samples and ds
  int num_samples = ceil(trajectory_length / ds_desired);
  Scalar ds = trajectory_length / ((Scalar)(num_samples - 1));
  // INFO_VAR(num_samples)
  // INFO_VAR(ds)

  DroneState beginstate = trajectory.state_in_time(0);
  samples.push_back(beginstate.p);

  Scalar tr_start_time = 0;
  Scalar tr_end_time = ttime;


  Scalar lower_time = 0;
  Scalar upper_time = 0;
  Scalar dist_first = 0;
  Scalar time_middle = 0;
  Scalar t_current = 0;
  Scalar t_from_start = 0;
  Scalar accumulated_s = 0;
  for (size_t i = 1; i < num_samples; i++) {
    const Scalar s = i * ds;
    Scalar target_ds = s - accumulated_s;

    // INFO("i " << i << " s " << s << " target_ds " << target_ds)
    lower_time = t_current;
    upper_time = tr_end_time;

    do {
      time_middle = (lower_time + upper_time) / 2.0;
      dist_first = trajectory.get_length(t_current, time_middle);

      if (dist_first > target_ds) {
        // soulution between lower_time and time_middle
        // INFO("bellow")
        upper_time = time_middle;
      } else {
        // soulution between time_middle and upper_time
        // INFO("above")
        lower_time = time_middle;
      }
    } while (fabs(dist_first - target_ds) > 0.001);

    /*
        // add also samples of acc switch
        std::vector<Scalar> acc_switch_times;
        for (size_t axi = 0; axi < 3; axi++) {
          const Scalar t1 = trajectory.get_axis_switch_time(axi);
          // const Scalar t1 = trajectory.get_axis(axi).t1;
          if (t1 > t_current && t1 < time_middle) {
            acc_switch_times.push_back(t1);
          }
        }
        std::sort(acc_switch_times.begin(), acc_switch_times.end());
    */

    // converged
    t_current = time_middle;
    accumulated_s += dist_first;

    // INFO("converged to time " << t_current);

    if (t_current < 0 || t_current > trajectory.timemin()) {
      INFO("bad time");
      exit(1);
    }


    DroneState dronestate = trajectory.state_in_time(t_current);

    const Scalar dist_to_last = (samples.back() - dronestate.p).norm();
    if (fabs(dist_to_last - ds) > ds) {
      INFO("bad distance change in sampleTrajectory " << dist_to_last)
      INFO("sample " << i << " out of " << num_samples)
      INFO_VAR(samples[0].size());
      INFO_VAR(target_ds);
      INFO_VAR(t_current)
      INFO_VAR(ttime)
      INFO_VAR(trajectory)
      INFO_VAR(dist_first)
      INFO_VAR(target_ds)
      INFO_VAR(fabs(dist_first - target_ds))
      exit(1);
    }
    samples.push_back(dronestate.p);
  }

  // INFO("sample trajectory end")
  return samples;
}


void VelocitySearchGraph::saveSamplesToFile(std::string filename,
                                            std::vector<QuadState> samples) {
  std::ofstream myfile;
  myfile.open(filename.c_str());
  if (myfile.is_open()) {
    if (samples.size() > 0 && samples[0].size() > 0) {
      myfile << to_string_raw_header(samples[0]) << std::endl;
    }

    const int s1 = samples.size();
    for (int var1 = 0; var1 < s1; ++var1) {
      const DroneState& data = samples[var1];
      myfile << to_string_raw(data);
      myfile << std::endl;
    }

    myfile.close();
  }
}

}  // namespace agi