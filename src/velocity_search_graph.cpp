
#include "velocity_search_graph.hpp"

#include <cfloat>
#include <fstream>
#include <iomanip>

#include "gravity.hpp"
#include "timer.hpp"

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

// Given waypoints that form a trajectory, finds velocities and angles for each
// waypoint such that the cost of trajectory traversal is optimal.

// Cone refocusing is used for sampling.
MultiWaypointTrajectory VelocitySearchGraph::find_velocities_in_positions(
  const std::vector<Vector<3>>& gates_waypoints,
  const Vector<3>& start_velocity, const Vector<3>& end_velocity,
  const std::vector<Scalar>& gates_yaw_deg,
  const std::vector<Scalar>& gates_pitch_deg,
  const std::vector<Scalar>& gates_vel_norms, const bool end_free,
  const bool use_gd) {


  const int gates_size = gates_waypoints.size();

  //  const std::vector<Scalar>& gates_yaw_deg,
  //  const std::vector<Scalar>& gates_pitch_deg,
  //  const std::vector<Scalar>& gates_vel_norms,
  //  const bool end_free


  // Scalar single_axis = sqrt(pow_max_acc_2 / 3.0);

//  for (auto x : gates_yaw_deg) {
//    std::cout << " HERE WP -> " << x << std::endl;
//  }
  // [Question#] What does this scalar represent? How is 'single_axis' max
  // acceleration different from 'mac_acc_'?
  Scalar pow_max_acc_2 = max_acc_ * max_acc_;
  // solution to a_s^2 + a_s^2 + (a_s+g)^2 = thrust^2
  Scalar single_axis =
    (-2 * G + sqrt(4 * G * G - 12 * (G * G - pow_max_acc_2))) / 6.0;

  Vector<3> max_acc_per_axis = Vector<3>::Constant(max_acc_);
//  Vector<3> max_acc_per_axis = Vector<3>::Constant(1.06);
  std::cout << "max_acc_ " << max_acc_ << std::endl;
//  std::cout << "single_axis " << single_axis << std::endl;

  if (gates_waypoints.size() != gates_yaw_deg.size() ||
      gates_yaw_deg.size() != gates_pitch_deg.size() ||
      gates_pitch_deg.size() != gates_vel_norms.size()) {
    std::cout << "wrong size of gates/yaw/pitch/norm vectors" << std::endl;
    return MultiWaypointTrajectory();
  }
  std::vector<Vector<3>> found_gates_speeds;
  // [Question#] What does found_gates_times[i] represent?
  std::vector<Scalar> found_gates_times;
  // Allow optimizing the end but omit the start. Therefore, "gates_size - 1"
  found_gates_speeds.resize(gates_size - 1);
  found_gates_times.resize(gates_size - 1);

  std::vector<std::vector<Eigen::Vector2f>> gates_yaw_pitch_size_ranges;
  gates_yaw_pitch_size_ranges.resize(
    gates_size - 1, {Eigen::Vector2f(-max_yaw_pitch_ang_, max_yaw_pitch_ang_),
                     Eigen::Vector2f(0, 0),
                     Eigen::Vector2f(min_velocity_size_, max_velocity_size_)});


  // [Question#] What is vel_size?
  // Samples have [gateid][sample_id] =
  //          (velocity_vector, (yaw, pitch,vel_size), (yaw_coneindex, pitch_coneindex, vel_size_coneindex))
  //
  std::vector<std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>>
    gate_velocity_samples;
  gate_velocity_samples.resize(gates_size);
  // For 'start' location, we don't create samples
  gate_velocity_samples[0].push_back(
    {start_velocity, Vector<3>::Zero(), Vector<3>::Zero()});


  // [Question#] What does time from the start mean?
  // shortest time structure for samples, add also the start, pair ==
  // {id_of_sample_in_gate_before,time from the start}
  std::vector<std::vector<std::pair<int, Scalar>>> shortest_samples_times;
  shortest_samples_times.resize(gate_velocity_samples.size());
  // For 'start' location, it takes 0 seconds to reach.
  // HOW TO INTERPRET --> To reach gate 0, it takes 0 seconds if we go from sample id '-1' from previous gate('-1')
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
  std::cout << "num_vel_samples_size -> " << num_vel_samples_size << std::endl;
  std::cout << "vel_samples_size_step -> " << vel_samples_size_step << std::endl;
  // We create sample for all gates except start and end
  std::cout << "gate_velocity_samples size -> " << gate_velocity_samples.size() << std::endl;
  std::cout << "shortest_samples_times size -> " << shortest_samples_times.size() << std::endl;

  // Amount of gates to create samples for. Initially, it's all gates in the middle
  // (all except start an end). Although, if 'end_free' = True, we also create samples
  // for end location.
  int sample_num_gates = gates_size - 2;
  if (end_free) {
    // If end is "free to sample", we have one more gate to sample velocities for.
    sample_num_gates += 1;
  } else {
    // For 'end' location, we don't create samples
    gate_velocity_samples[gates_size - 1].push_back(
      {end_velocity, Vector<3>::Zero(), Vector<3>::Zero()});
    // It takes infinite amount time to reach to the end gate.
    shortest_samples_times[gate_velocity_samples.size() - 1].push_back(
      {-1, DBL_MAX});  // add end
  }

  // change the velocity norm range based on current velocity in the
  // gates_vel_norms
  for (size_t gid = 1; gid <= sample_num_gates; gid++) {
    // minimal distance from norm boundaries to desired velocity
    Eigen::Vector2f& min_max_norm = gates_yaw_pitch_size_ranges[gid - 1][2];
    const Scalar gate_vel_norm = gates_vel_norms[gid];
    // std::cout << "gate_vel_norm " << gate_vel_norm << std::endl;
    Scalar min_dist_desired_vel = std::min(gate_vel_norm - min_max_norm(0),
                                           min_max_norm(1) - gate_vel_norm);
    min_dist_desired_vel = std::max(min_dist_desired_vel, 2.0);
    min_max_norm(0) = gate_vel_norm - min_dist_desired_vel;
    min_max_norm(1) = gate_vel_norm + min_dist_desired_vel;
  }

  Timer timer;
  timer.tic();
  Timer timer_samples;
  Timer timer_calc_pm;

  const int num_iter_improve = 10;

  auto backup_found_gates_speeds = found_gates_speeds;
  auto backup_gate_velocity_samples = gate_velocity_samples;
  auto backup_shortest_samples_times = shortest_samples_times;

  // Shortest time of the whole trajectory.
  Scalar last_shortest_time = MAX_SCALAR;
  // iteration over improvements of velocities == focusing velocity cones
  for (size_t iterimp = 0; iterimp < num_iter_improve; iterimp++) {

    backup_found_gates_speeds = found_gates_speeds;
    backup_gate_velocity_samples = gate_velocity_samples;
    backup_shortest_samples_times = shortest_samples_times;

    // COPY BACK TO => found_gates_speeds, gate_velocity_samples, shortest_samples_times

    // loop gates 1....sample_num_gates (0 is start and gatesize-1 is end)
    std::cout << "Improvement -> " << iterimp << std::endl;
    timer_samples.tic();
    // here creating only the samples
    for (size_t gid = 1; gid <= sample_num_gates; gid++) {
      const int num_samples_gate =
        num_vel_samples_yaw_pitch_ang * num_vel_samples_yaw_pitch_ang *
        num_vel_samples_size;  // zero size does not need yaw-pitch samples

      // get current min-max angles(relative to gate yaw) and velocity sizes
      const Eigen::Vector2f min_max_yaw =
        gates_yaw_pitch_size_ranges[gid - 1][0];
      const Eigen::Vector2f min_max_pitch =
        gates_yaw_pitch_size_ranges[gid - 1][1];
      const Eigen::Vector2f min_max_norm =
        gates_yaw_pitch_size_ranges[gid - 1][2];

      // angle step sizes
      const Scalar gate_vel_samples_size_step =
        (min_max_norm(1) - min_max_norm(0)) /
        ((Scalar)(num_vel_samples_size - 1));
      const Scalar gate_vel_samples_yaw_step =
        (min_max_yaw(1) - min_max_yaw(0)) /
        ((Scalar)(num_vel_samples_yaw_pitch_ang - 1));
      const Scalar gate_vel_samples_pitch_step =
        (min_max_pitch(1) - min_max_pitch(0)) /
        ((Scalar)(num_vel_samples_yaw_pitch_ang - 1));

      const Scalar gate_yaw = gates_yaw_deg[gid];
      const Scalar gate_pitch = gates_pitch_deg[gid];
      const Scalar gate_vel_norm = gates_vel_norms[gid];

      // For gate 'gid' there is going to be 'num_samples_gate' samples.
      gate_velocity_samples[gid].resize(num_samples_gate);
      // For gate 'gid' there is going to be 'num_samples_gate' samples, and it would
      // Take infinite amonut of time to reach that gate, initially.
      shortest_samples_times[gid].resize(num_samples_gate, {-1, DBL_MAX});
      const std::pair<int, Scalar> not_reached = std::make_pair(-1, DBL_MAX);
      std::fill(shortest_samples_times[gid].begin(),
                shortest_samples_times[gid].end(), not_reached);

//            for (auto x : shortest_samples_times[gid]) {
//              std::cout << "[0] -> " << x.first << ", " << "[1] -> " << x.second << std::endl;
//            }

      // save samples values to the gate_velocity_samples
      int i = 0;
      for (size_t z_s = 0; z_s < num_vel_samples_size; z_s++) {
        const Scalar vel_size =
          min_max_norm(0) + z_s * gate_vel_samples_size_step;
        for (size_t x_s = 0; x_s < num_vel_samples_yaw_pitch_ang; x_s++) {
          const Scalar yaw = min_max_yaw(0) + x_s * gate_vel_samples_yaw_step;
          for (size_t y_s = 0; y_s < num_vel_samples_yaw_pitch_ang; y_s++) {
            const Scalar pitch =
              min_max_pitch(0) + y_s * gate_vel_samples_pitch_step;
            Matrix<3, 3> m;
            m = Eigen::AngleAxisd((gate_pitch + pitch) * M_PI / 180.0,
                                  Vector<3>::UnitY()) *
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
    timer_samples.toc();


    // Now we have the samples, search the shortest time path through the graph.
    // [Question#] We take two consecutive points and find the shortest distance between those?


    // ---- FIND OPTIMAL VELOCITIES AND TIME OF THE TRAJECTORY.

    // Here we start finding velocity in points
    // For each point in trajectory.
    for (size_t gid_to = 1; gid_to < gates_size; gid_to++) {
//       std::cout << "gid_to " << gid_to << std::endl;
      // get locally the parts
      const std::vector<std::pair<int, Scalar>>& from_shortest_samples_dists =
        shortest_samples_times[gid_to - 1];
      std::vector<std::pair<int, Scalar>>& to_shortest_samples_dists =
        shortest_samples_times[gid_to];
//      for (auto x : to_shortest_samples_dists) {
//        std::cout << "gid_to " << gid_to << ", Pair" << "<" << x.first << ", " << x.second << "> \n";
//      }
      const std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>&
        from_samples = gate_velocity_samples[gid_to - 1];
      const std::vector<std::tuple<Vector<3>, Vector<3>, Vector<3>>>&
        to_samples = gate_velocity_samples[gid_to];

      // get the positions of the gates/start/end
      const Vector<3>& from_p = gates_waypoints[gid_to - 1];
      const Vector<3>& to_p = gates_waypoints[gid_to];

      Scalar min_dist_to_gate = DBL_MAX;

      // loop for indexes of the samples in the gate (gid_to-1)
      for (int idx_from = 0; idx_from < from_shortest_samples_dists.size(); idx_from++) {
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
          PointMassTrajectory3D tr_max_acc(from_state, to_state, max_acc_per_axis, true);

          if (tr_max_acc.exists()) {
//            std::cout << "from " << from_state.p.transpose() << " "
//                      << from_state.v.transpose() << std::endl;
//            std::cout << "to " << to_state.p.transpose() << " "
//                      << to_state.v.transpose() << std::endl;
//            std::cout << tr_max_acc << std::endl;
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

//      for (auto x : shortest_samples_times[2]) {
//        std::cout << "gid_to " << 2 << ", Pair" << "<" << x.first << ", " << x.second << "> \n";
//      }

      if (min_dist_to_gate == DBL_MAX) {
        std::cout << "there is no connection to gate " << gid_to
                  << " with position " << gates_waypoints[gid_to].transpose()
                  << std::endl;
//        std::cout << from_p.transpose() << std::endl;
//        std::cout << to_p.transpose() << std::endl;
        return MultiWaypointTrajectory();
      }
    }
    // ---- FIND OPTIMAL VELOCITIES AND TIME OF THE TRAJECTORY. ---- END

    // Find the time of the shortest-time trajectory by looking at samples of
    // the last node
//    std::cout << "Heyaaa ----------------" << std::endl;
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
//    std::cout << "end_best_idx -> " << end_best_idx << std::endl;

    std::cout << "last_shortest_time blah -> " << last_shortest_time << std::endl;
    std::cout << "shortest_time blah -> " << shortest_time << std::endl << std::endl;
    // We continue this iterative
    // process until the ratio between previous minimum time and
    // new minimum time is less than . In our case, this  has
    // been chosen such that if the improvement is less than 1% in
    // minimum time between iteration, the algorithm has converged
    // and the iteration ends.
    if (shortest_time > 0.99 * last_shortest_time) {
      std::cout << "We break" << std::endl;
      break;
    }
    last_shortest_time = shortest_time;

    // [Question#] What are we trying to optimize here?
    // additionally also optimize the ranges of the velocity cone
    // trace back the trajectory from the shortest final sample

    // Knowing the id of the sample of the last location that gives us the best
    // time, we can traceback the optimal velocities at all previous locations.
    int prev_sample_idx = end_best_idx;
    for (int g_id = endi; g_id > 0; g_id--) {
      found_gates_speeds[g_id - 1] =
        std::get<0>(gate_velocity_samples[g_id][prev_sample_idx]);
//      std::cout << "idx " << g_id << " speed "
//                << found_gates_speeds[g_id - 1].transpose() << std::endl;
      Vector<3> indexes =
        std::get<2>(gate_velocity_samples[g_id][prev_sample_idx]);
      for (size_t i = 0; i < 3; i++) {
        const Eigen::Vector2f current_range =
          gates_yaw_pitch_size_ranges[g_id - 1][i];

        const Scalar range_size_half =
          (current_range(1) - current_range(0)) / 2.0;
        if (indexes(i) == 0) {
          // move the range of minmax down
          gates_yaw_pitch_size_ranges[g_id - 1][i] =
            Eigen::Vector2f(current_range(0) - range_size_half,
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
            Eigen::Vector2f(current_range(0) + range_size_half,
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
            Eigen::Vector2f(current_range(0) + range_size_half / 2.0,
                            current_range(1) - range_size_half / 2.0);
        }
//        std::cout << "g " << (g_id - 1) << " ax " << i << "range"
//                  << gates_yaw_pitch_size_ranges[g_id - 1][i].transpose()
//                  << std::endl;
      }

      // traceback previous sample index here
      prev_sample_idx = shortest_samples_times[g_id][prev_sample_idx].first;
    }
  }
  std::cout << "Exiting the iter improvement loop.." << std::endl;

  timer.toc();

  // COPY BACK TO => found_gates_speeds, gate_velocity_samples, shortest_samples_times
  found_gates_speeds = backup_found_gates_speeds;
  gate_velocity_samples = backup_gate_velocity_samples;
  shortest_samples_times = backup_shortest_samples_times;


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

   std::cout << "last_shortest_time -> " << last_shortest_time << std::endl;

  int prev_sample_idx = end_best_idx;
  for (int g_id = endi; g_id > 0; g_id--) {
    found_gates_speeds[g_id - 1] =
      std::get<0>(gate_velocity_samples[g_id][prev_sample_idx]);
    found_gates_times[g_id - 1] =
      shortest_samples_times[g_id][prev_sample_idx].second;

    prev_sample_idx = shortest_samples_times[g_id][prev_sample_idx].first;
  }


  // Populate trajectories between all points.
  MultiWaypointTrajectory trajectories;
  // Total time of all trajectories.
  Scalar time_sum = 0;

  // ------ Add trajectory "start" -> "location_1"
  QuadState from_start_state;
  from_start_state.p = gates_waypoints[0];
  from_start_state.v = start_velocity;
  QuadState to_state;
  to_state.p = gates_waypoints[1];
  to_state.v = found_gates_speeds[0];
  PointMassTrajectory3D tr_max_acc_start;
  tr_max_acc_start =
    PointMassTrajectory3D(from_start_state, to_state, max_acc_per_axis, true);
  trajectories.push_back(tr_max_acc_start);
  time_sum += tr_max_acc_start.time();
  if (!tr_max_acc_start.exists()) {
    std::cout << "Trajectory does not exist from start location." << std::endl;
    return MultiWaypointTrajectory();
  }
  // ------ Add trajectory "start" -> "location_1" ---- END

  // ------ Add all other trajectories
  for (size_t i = 2; i < gates_waypoints.size(); i++) {
    QuadState from_state_b;
    from_state_b.p = gates_waypoints[i - 1];
    from_state_b.v = found_gates_speeds[i - 2];
    QuadState to_state_b;
    to_state_b.p = gates_waypoints[i];
    to_state_b.v = found_gates_speeds[i - 1];

    PointMassTrajectory3D tr_max_between;
    tr_max_between =
      PointMassTrajectory3D(from_state_b, to_state_b, max_acc_per_axis, true);
    trajectories.push_back(tr_max_between);
    time_sum += tr_max_between.time();

    if (!tr_max_between.exists()) {
      // This trajectory does not exist.
      std::cout << "found pmm trajectory does not exist, gate " << i
                << std::endl;
      std::cout << "equalized" << std::endl;
      std::cout << tr_max_between << std::endl;
      return MultiWaypointTrajectory();
    }
    if (fabs(tr_max_between.time_min() - tr_max_between.time()) >=
        PRECISION_PMM_VALUES) {
      std::cout << tr_max_between << std::endl;
      std::cout << "not equal time " << std::endl;
//      exit(1);
    }
  }
  // ------ Add all other trajectories --- END
  if (fabs(time_sum - shortest_time) > 0.0001) {
    std::cout << "error calculating pmm trajectory" << std::endl;
    std::cout << "time_sum " << time_sum << std::endl;
    std::cout << "shortest_time" << shortest_time << std::endl;
    return MultiWaypointTrajectory();
  }
  Scalar actual_time = 0;
  for (auto tr: trajectories) {
    actual_time += tr.time();
  }
//  std::cout << "End of trajectory calculation. Return trajectories found.." << std::endl;
//  std::cout << "Actual time -> " << actual_time << std::endl;
  std::cout << "Actual time -> " << time_sum << std::endl;
//  for (auto t : found_gates_times) {
//    std::cout << "time -> " << t << std::endl;
//  }
  return trajectories;
}

void VelocitySearchGraph::setYawPitchConeSidesAndPecisiton(
  const Scalar max_yaw_pitch_ang, const Scalar precision_yaw_pitch_ang) {
  max_yaw_pitch_ang_ = max_yaw_pitch_ang;
  precision_yaw_pitch_ang_ = precision_yaw_pitch_ang;
}

void VelocitySearchGraph::saveTrajectoryEquidistant(
  const MultiWaypointTrajectory& trajectories, std::string filename,
  const Scalar desired_ds) {
  // std::cout << "saveTrajectoryEquidistant " << desired_ds << std::endl;
  std::vector<QuadState> samples =
    getTrajectoryEquidistantStates(trajectories, desired_ds);
  // std::cout << "have samples " << std::endl;
  saveSamplesToFile(filename, samples);
}

std::vector<QuadState> VelocitySearchGraph::getTrajectoryEquidistantStates(
  const MultiWaypointTrajectory& trajectories, const Scalar desired_ds) {
  std::vector<QuadState> samples;
  Scalar trajectory_length = 0;
//  std::cout << "here1" << std::endl;
  // get the trajectories length
  for (size_t i = 0; i < trajectories.size(); i++) {
    const Scalar ttime = trajectories[i].time_min();
    const Scalar traj_length =
      trajectories[i].get_length_between_times(0, ttime);
    trajectory_length += traj_length;
  }
//  std::cout << "here" << std::endl;
  const int num_samples = ceil(trajectory_length / desired_ds);
  const Scalar ds = trajectory_length / ((Scalar)(num_samples - 1));

  QuadState beginstate = trajectories[0].state_in_time(0);
  beginstate.t = 0;
  samples.push_back(beginstate);

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
  for (size_t i = 1; i < num_samples; i++) {
    const Scalar s = i * ds;
    // std::cout << "num_samples " << num_samples << " i " << i << std::endl;
    Scalar target_ds = s - accumulated_s;
    // std::cout << "target_ds " << target_ds << std::endl;
    Scalar dist_to_end =
      trajectories[tr_id].get_length_between_times(t_current, tr_end_time);
    // std::cout << "dist_to_end " << dist_to_end << std::endl;

    while (dist_to_end < target_ds) {
      tr_id++;
      if (tr_id >= trajectories.size()) {
        break;
      }
      t_current = 0;
      t_from_start += tr_end_time;
      tr_end_time = trajectories[tr_id].time_min();
      target_ds -= dist_to_end;
      accumulated_s += dist_to_end;
      dist_to_end =
        trajectories[tr_id].get_length_between_times(t_current, tr_end_time);
    }

    if (tr_id >= trajectories.size()) {
      break;
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
    // std::cout << "dist_first " << dist_first << std::endl;
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
      return std::vector<QuadState>();
    }

    QuadState dronestate = trajectories[tr_id].state_in_time(t_current);
    dronestate.t = t_from_start + t_current;
    /*
    const Scalar dist_to_last = (samples.back().p - dronestate.p).norm();
    if (fabs(dist_to_last - ds) > 0.10 * ds) {
      std::cout << "bad distance change in save_track_trajectory_equidistant "
                << dist_to_last << std::endl;
      std::cout << "ds " << ds << std::endl;
      std::cout << "sample " << i << " out of " << num_samples << std::endl;
      return std::vector<QuadState>();
    }
    */
    samples.push_back(dronestate);
  }

  return samples;
}

void VelocitySearchGraph::saveTrajectoryEquitemporal(
  const MultiWaypointTrajectory& trajectories, std::string filename,
  const Scalar dt_desired) {
  std::vector<QuadState> samples;

  Scalar trajectories_time = 0;
  for (size_t i = 0; i < trajectories.size(); i++)
    trajectories_time += trajectories[i].time();

  const int num_samples = ceil(trajectories_time / dt_desired);

  samples.resize(num_samples + 1);
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
    samples[i] = ds;
  }

  saveSamplesToFile(filename, samples);
}

std::string VelocitySearchGraph::stateToStringHeader(const QuadState& s) {
  std::stringstream o;
  o << std::setprecision(4)
    << "t,p_x,p_y,p_z,q_w,q_x,q_y,q_z,v_x,v_y,v_z,w_x,w_y,w_z,a_lin_x,a_lin_y,"
       "a_lin_z";
  return o.str();
}

std::string VelocitySearchGraph::stateToString(const QuadState& s) {
  std::stringstream o;
  o << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << s.t
    << "," << s.p(0) << "," << s.p(1) << "," << s.p(2) << "," << s.qx(0) << ","
    << s.qx(1) << "," << s.qx(2) << "," << s.qx(3) << "," << s.v(0) << ","
    << s.v(1) << "," << s.v(2) << "," << s.w(0) << "," << s.w(1) << ","
    << s.w(2) << "," << s.a(0) << "," << s.a(1) << "," << s.a(2);
  return o.str();
}

void VelocitySearchGraph::saveSamplesToFile(std::string filename,
                                            std::vector<QuadState> samples) {
  std::ofstream myfile;
  myfile.open(filename.c_str());
  if (myfile.is_open()) {
    if (samples.size() > 0 && samples[0].size() > 0) {
      myfile << stateToStringHeader(samples[0]) << std::endl;
    }

    const int s1 = samples.size();
    for (int var1 = 0; var1 < s1; ++var1) {
      const QuadState& data = samples[var1];
      myfile << stateToString(data);
      myfile << std::endl;
    }

    myfile.close();
  }
}

void operator>>(const YAML::Node& node, Scalar& value) {
  assert(node.IsScalar());
  value = node.as<Scalar>();
}

void operator>>(const YAML::Node& node, agi::Vector<3>& value) {
  for (std::size_t i = 0; i < node.size(); i++) {
    value(i) = node[i].as<Scalar>();
  }
}

}  // namespace agi