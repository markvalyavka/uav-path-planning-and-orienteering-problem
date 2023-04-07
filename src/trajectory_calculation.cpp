
#include "trajectory_calculation.hpp"


namespace agi {


std::tuple<MultiWaypointTrajectory, Scalar> calculate_trajectory_cost_and_optimal_velocities(std::vector<int> &scheduled_locations_idx,
                                                                                             EnvConfig& env_params,
                                                                                             bool sample_start_velocity) {

  travel_cost_map& travel_costs = env_params.precalculated_costs;
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>>& velocity_sample_tuples = env_params.velocity_samples_tuples;
  Vector<3> &start_velocity = env_params.start_velocity;
  Vector<3> &max_acc_per_axis = env_params.max_acc_per_axis;
  std::vector<Vector<3>> location_positions{};
  for (int idx : scheduled_locations_idx) {
    location_positions.push_back(env_params.location_positions[idx]);
    //    std::cout << idx << std::endl;
  }

  Scalar shortest_time = DBL_MAX;
  int number_of_locations = location_positions.size();
  int number_of_velocity_samples = velocity_sample_tuples.size();
  std::tuple<Scalar, Scalar> start_vel_mag_ang = to_velocity_norm_and_angle(start_velocity);

  // {sample_id_in_location_before, total_time_from_the_start}, idx - location_idx
  std::vector<std::vector<std::pair<int, Scalar>>> shortest_samples_times;
  shortest_samples_times.resize(number_of_locations);


  const std::pair<int, Scalar> not_reached = std::make_pair(-1, DBL_MAX);
  for (int loc_id = 1; loc_id < number_of_locations; loc_id++) {
    // Every sample at `loc_id` is unreachable initially.
    shortest_samples_times[loc_id].resize(number_of_velocity_samples, {-1, DBL_MAX});
    std::fill(shortest_samples_times[loc_id].begin(),
              shortest_samples_times[loc_id].end(), not_reached);
  }
  // For 'start' location, it takes 0 seconds to reach.
  if (sample_start_velocity) {
    const std::pair<int, Scalar> reached = std::make_pair(-1, 0);
    shortest_samples_times[0].resize(number_of_velocity_samples, {-1, 0});
    std::fill(shortest_samples_times[0].begin(),
              shortest_samples_times[0].end(), reached);
  } else {
    shortest_samples_times[0].push_back({-1, 0});
  }


  // Samples have [loc_id][sample_id][{velocity, norm, heading_angle)}]
  // number of `sample_ids` == V * H, essentially.
  std::vector<std::vector<std::tuple<Vector<3>, Scalar, Scalar>>> gate_velocity_samples;
  gate_velocity_samples.resize(number_of_locations);
  for (int loc_id = 1; loc_id < number_of_locations; loc_id++) {
    gate_velocity_samples[loc_id] = velocity_sample_tuples;
  }
  if (sample_start_velocity) {
    gate_velocity_samples[0] = velocity_sample_tuples;
  } else {
    gate_velocity_samples[0].push_back({start_velocity, 0, 0});
  }


  // ---- FIND OPTIMAL VELOCITIES AND TIME OF THE TRAJECTORY.

  for (size_t loc_to = 1; loc_to < number_of_locations; loc_to++) {
    //    std::cout << "loc_to " << loc_to << std::endl;
    const std::vector<std::pair<int, Scalar>>& from_shortest_samples_dists =
      shortest_samples_times[loc_to - 1];
    std::vector<std::pair<int, Scalar>>& to_shortest_samples_dists =
      shortest_samples_times[loc_to];
    const std::vector<std::tuple<Vector<3>, Scalar, Scalar>>&
      from_samples = gate_velocity_samples[loc_to - 1];
    const std::vector<std::tuple<Vector<3>, Scalar, Scalar>>&
      to_samples = gate_velocity_samples[loc_to];

    // get the positions of the gates/start/end
    const Vector<3>& from_p = location_positions[loc_to - 1];
    const Vector<3>& to_p = location_positions[loc_to];

    Scalar min_dist_to_gate = DBL_MAX;

    // loop for indexes of the samples in the gate (gid_to-1)
    for (int sample_idx_from = 0; sample_idx_from < from_shortest_samples_dists.size(); sample_idx_from++) {
      const Scalar& time_from = from_shortest_samples_dists[sample_idx_from].second;
      const Vector<3>& from_v = std::get<0>(from_samples[sample_idx_from]);
      if (time_from == DBL_MAX) {
        // No path reaching that sample.
        continue;
      }
      // loop for indexes of the samples in the gate gid_to
      for (int sample_idx_to = 0; sample_idx_to < to_shortest_samples_dists.size(); sample_idx_to++) {
        std::pair<int, Scalar>& shortest_to = to_shortest_samples_dists[sample_idx_to];

        Scalar from_norm = std::get<1>(from_samples[sample_idx_from]);
        Scalar from_angle = std::get<2>(from_samples[sample_idx_from]);
        Scalar to_norm = std::get<1>(to_samples[sample_idx_to]);
        Scalar to_angle = std::get<2>(to_samples[sample_idx_to]);

        // If location is START and the velocity is fixed (not sampled).
        if (!sample_start_velocity && loc_to-1 == 0) {
          from_norm = std::get<0>(start_vel_mag_ang);
          from_angle = std::get<1>(start_vel_mag_ang);
          //          std::cout << "blah" << std::endl;
        }

        Scalar time_between = travel_costs[scheduled_locations_idx[loc_to-1]][scheduled_locations_idx[loc_to]][from_norm][to_norm][from_angle][to_angle];

        const Scalar time_tot = time_between + time_from;
        if (time_tot < min_dist_to_gate) {
          min_dist_to_gate = time_tot;
        }
        if (time_tot < shortest_to.second) {
          shortest_to.first = sample_idx_from;
          shortest_to.second = time_tot;
        }
      }
    }
    if (min_dist_to_gate == DBL_MAX) {
      std::cout << "Location " << loc_to
                << "[" << location_positions[loc_to].transpose() << "]"
                << " is unreachable." << std::endl;
      return {MultiWaypointTrajectory(), DBL_MAX};
    }
  }
  // ---- FIND OPTIMAL VELOCITIES AND TIME OF THE TRAJECTORY. ---- END

  int end_best_idx = -1;
  const int endi = shortest_samples_times.size() - 1;
  for (size_t i = 0; i < shortest_samples_times[endi].size(); i++) {
    const Scalar time = shortest_samples_times[endi][i].second;
    if (time < shortest_time) {
      shortest_time = time;
      end_best_idx = i;
    }
  }

  //  std::cout << "end_best_idx -> " << end_best_idx << std::endl;
  //    std::cout << "shortest_time " << shortest_time << std::endl << std::endl;

  // Knowing the id of the sample of the last location that gives us the best
  // time, we can traceback the optimal velocities at all previous locations.
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>> found_optimal_gates_velocities;
  found_optimal_gates_velocities.resize(number_of_locations);
  int prev_sample_idx = end_best_idx;
  for (int loc_id = endi; loc_id >= 0; loc_id--) {
    found_optimal_gates_velocities[loc_id] = gate_velocity_samples[loc_id][prev_sample_idx];
    //    found_gates_times[loc_id] = shortest_samples_times[loc_id][prev_sample_idx].second;

    prev_sample_idx = shortest_samples_times[loc_id][prev_sample_idx].first;
  }
  if (!sample_start_velocity) {
    found_optimal_gates_velocities[0] = {start_velocity, std::get<0>(start_vel_mag_ang), std::get<1>(start_vel_mag_ang)};
  }


  // Populate trajectories between all points.
  MultiWaypointTrajectory trajectories;
  // Total time of all trajectories.
  Scalar time_sum = 0;

  // ------ Add trajectory "start" -> "location_1"
  //  QuadState from_start_state;
  //  from_start_state.p = location_positions[0];
  //  from_start_state.v = start_velocity;
  //  from_start_state.v_norm = std::get<0>(start_vel_mag_ang);
  //  from_start_state.v_angle = std::get<1>(start_vel_mag_ang);;
  //  QuadState to_state;
  //  to_state.p = location_positions[1];
  //  to_state.v = std::get<0>(found_optimal_gates_velocities[0]);
  //  to_state.v_norm = std::get<1>(found_optimal_gates_velocities[0]);
  //  to_state.v_angle = std::get<2>(found_optimal_gates_velocities[0]);
  //  PointMassTrajectory3D tr_max_acc_start;
  //  tr_max_acc_start =
  //    PointMassTrajectory3D(from_start_state, to_state, max_acc_per_axis, true);

  //  trajectories.push_back(tr_max_acc_start);
  //  time_sum += tr_max_acc_start.time();
  //  if (!tr_max_acc_start.exists()) {
  //    std::cout << "Trajectory does not exist from start location." << std::endl;
  //    return {MultiWaypointTrajectory(), DBL_MAX};
  //  }
  // ------ Add trajectory "start" -> "location_1" ---- END

  // ------ Add all other trajectories
  for (size_t i = 1; i < location_positions.size(); i++) {
    QuadState from_state_b;
    from_state_b.p = location_positions[i - 1];
    from_state_b.v = std::get<0>(found_optimal_gates_velocities[i - 1]);
    from_state_b.v_norm = std::get<1>(found_optimal_gates_velocities[i - 1]);
    from_state_b.v_angle = std::get<2>(found_optimal_gates_velocities[i - 1]);

    QuadState to_state_b;
    to_state_b.p = location_positions[i];
    to_state_b.v = std::get<0>(found_optimal_gates_velocities[i]);
    to_state_b.v_norm = std::get<1>(found_optimal_gates_velocities[i]);
    to_state_b.v_angle = std::get<2>(found_optimal_gates_velocities[i]);

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
      return {MultiWaypointTrajectory(), DBL_MAX};
    }
    if (fabs(tr_max_between.time_min() - tr_max_between.time()) >=
        PRECISION_PMM_VALUES) {
      //      std::cout << tr_max_between << std::endl;
      //      std::cout << fabs(tr_max_between.time_min() - tr_max_between.time()) << " >= " << PRECISION_PMM_VALUES << std::endl;
      //      std::cout << "Not equal time here" << std::endl;
    }
  }
  // ------ Add all other trajectories --- END

  if (fabs(time_sum - shortest_time) > 0.0001) {
    std::cout << "error calculating pmm trajectory" << std::endl;
    std::cout << "time_sum " << time_sum << std::endl;
    std::cout << "shortest_time" << shortest_time << std::endl;
    return {MultiWaypointTrajectory(), DBL_MAX};
  }
  return {trajectories, shortest_time};
}



}