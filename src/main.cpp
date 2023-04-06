#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <csignal>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <cstdlib>
#include <random>

#include "gravity.hpp"
#include "pmm_trajectory3d.hpp"
#include "three_acc.hpp"
#include "timer.hpp"
#include "tuples_hash.hpp"
#include "velocity_search_graph.hpp"
#include "env_config.hpp"
#include "yaml-cpp/yaml.h"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received. Exiting.. " << std::endl;
  exit(sig);
}

typedef std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>, std::vector<int>> constructed_trajectory;

int randint(int Min, int Max) {
  return std::rand() % (Max + 1 - Min) + Min;
}

MultiWaypointTrajectory test_pmm(std::vector<Vector<3>> gates_waypoints,
             std::vector<Scalar> gates_vel_norms,
             std::vector<Scalar> gates_yaw_deg,
             Vector<3> start_velocity) {
  const Scalar yaw_pitch_cone_angle_boundary = 40.0;
  // Velocity norm samples.
  const Scalar min_velocity_norm = 0;
  const Scalar min_velocity_norm_boundary = 0.0001;
  const Scalar max_velocity_norm = 2.111;
  // Distance between velocity norms
  const Scalar precision_velocity_norm = 0.1;
  const Scalar max_acc_norm = 1.06066;
  const bool end_free = true;
  Scalar max_yaw_pitch_ang = 20;
  Scalar precision_yaw_pitch_ang = 20;

  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  std::vector<Scalar> gates_pitch_deg(gates_yaw_deg.size(), 0);
  Vector<3> end_velocity(0, 0, 0);

  MultiWaypointTrajectory tr = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, end_velocity, gates_yaw_deg,
    gates_pitch_deg, gates_vel_norms, end_free, false);

  return tr;
}

  // point mass model
// Keep this code for code refocusing example
MultiWaypointTrajectory test_pmm() {

  // PARAMS:
  // 1. gates_yaw_deg
  // 2. gates_vel_norms
  // 3. gates_waypoints
  // 4. start_velocity


  // Register signal for killing
  std::signal(SIGINT, signal_callback);
  std::cout << "Testing PMM trajectories " << std::endl;

  std::string config_file = "config.yaml";
  YAML::Node config = YAML::LoadFile(config_file);

  // PARAMETERS INITIALIZATION

  // -+ Value of initial samples.
//  Scalar max_yaw_pitch_ang = 20.0;
  // Distance between initial samples.
//  Scalar precision_yaw_pitch_ang = 20.0;
  // Maximal values of the samples. (kinda leeway?)
  const Scalar yaw_pitch_cone_angle_boundary = 40.0;
  // Velocity norm samples.
  const Scalar min_velocity_norm = 0;
  const Scalar min_velocity_norm_boundary = 0.0001;
  const Scalar max_velocity_norm = 2.111;
  // Distance between velocity norms
  const Scalar precision_velocity_norm = 0.1;
  const Scalar max_acc_norm = 1.06066;
  const bool end_free = true;
  Scalar max_yaw_pitch_ang = 20;
  Scalar precision_yaw_pitch_ang = 20;

  // Construct a graph given all the parameters of the environment.
  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  Vector<3> start_velocity;
  Vector<3> start_position;
  Vector<3> end_velocity(0, 0, 0);
  Vector<3> end_position;
  config["start"]["velocity"] >> start_velocity;
  config["start"]["position"] >> start_position;
//  config["end"]["velocity"] >> end_velocity;
  config["end"]["position"] >> end_position;

  // Vector of 'heading angles'(positions) of all gates (including start and end).
  std::vector<Vector<3>> gates_waypoints;
  if (!parseArrayParam<Vector<3>>(config, "locations", gates_waypoints))
    std::cerr << "can't load param gates" << std::endl;
  gates_waypoints.insert(gates_waypoints.begin(), start_position);
  gates_waypoints.push_back(end_position);
  for (auto loc : gates_waypoints) {
    std::cout << loc.transpose() << std::endl;
  }
  std::vector<Scalar> gates_yaw_deg{90, 112, 112, 112, 67, 45, 22, 22, 180, 180, 22, 180, 180};
//  if (!parseArrayParam<Scalar>(config, "gates_orientations", gates_yaw_deg))
//    std::cerr << "can't load param gates_orientations" << std::endl;

  std::vector<Scalar> gates_pitch_deg(gates_yaw_deg.size(), 0);
//  if (!parseArrayParam<Scalar>(config, "gates_pitch_orientations",
//                               gates_pitch_deg))
//    std::cerr << "can't load param gates_pitch_orientations" << std::endl;
  std::vector<Scalar> gates_vel_norms{1.06066, 1.06066, 2.12132, 2.12132, 2.12132, 2.12132, 1.06066, 0, 1.06066, 2.12132, 0, 1.06066, 2.12132};
//  std::vector<Scalar> gates_vel_norms{};
//  gates_vel_norms.resize(gates_pitch_deg.size(),
//                         (max_velocity_norm + min_velocity_norm) / 2.0);



//  std::cout << "start_velocity " << start_velocity.transpose() << std::endl;
//  std::cout << "gates_waypoints.size() " << gates_waypoints.size() << std::endl;
//  std::cout << "gates_yaw_deg.size() " << gates_yaw_deg.size() << std::endl;
//  std::cout << "gates_pitch_deg.size() " << gates_pitch_deg.size() << std::endl;
//  std::cout << "gates_vel_norms.size() " << gates_vel_norms.size() << std::endl;



  Timer find_vel("find vel");
  Scalar sum_times = 0;
  find_vel.tic();


  MultiWaypointTrajectory tr;
  // Finds velocities, (and pitch and yaw?) in each gate that optimize the trajectory
  // to traverse all the gates. Plus, it gives information about the time to reach
  // from gateA to gateB.
  tr = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, end_velocity, gates_yaw_deg,
    gates_pitch_deg, gates_vel_norms, end_free, false);
  return tr;
//  for (auto t: tr) {
//    std::cout << t.inp_from_v_norm << std::endl;
//  }

  for (size_t i = 0; i < tr.size(); i++) {
    const Vector<3> vel = tr[i].get_end_state().v;
    const Vector<3> vel_norm = vel.normalized();
//    std::cout << "vel -> " << vel.transpose() << std::endl;
//    std::cout << "vel_norm -> " << vel.norm() << std::endl;
    const Scalar pitch = asin(-vel_norm(2)) * 180 / M_PI;
    const Scalar yaw = atan2(vel_norm(1), vel_norm(0)) * 180 / M_PI;
    gates_yaw_deg[i + 1] = yaw;
    gates_pitch_deg[i + 1] = pitch;
    gates_vel_norms[i + 1] = vel.norm();
  }

  for (size_t i = 0; i < tr.size(); i++) {
    Vector<3> vel_norm = tr[i].get_end_state().v.normalized();
    if (i == 0) {
      Vector<3> vel_norm = tr[i].get_start_state().v.normalized();
    }
    Scalar pitch = asin(-vel_norm(2)) * 180 / M_PI;
    Scalar yaw = atan2(vel_norm(1), vel_norm(0)) * 180 / M_PI;
//    std::cout << i + 1 << " pos " << tr[i].get_end_state().p.transpose()
//              << " vel " << tr[i].get_end_state().v.transpose() << " acc "
//              << tr[i].get_end_state().a.transpose() << " thrust norm "
//              << (tr[i].get_end_state().a - GVEC).norm() << " exists "
//              << tr[i].exists() << " yaw " << yaw << " pitch " << pitch
//              << std::endl;

    sum_times += tr[i].time();
//    if (tr[i].time() - tr[i].time_min() > PRECISION_PMM_VALUES) {
//      std::cout << "bad time!!!!!!" << std::endl;
//    }
  }
  std::cout <<  "--------------------" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(tr, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;

//  return 0;
}



std::tuple<MultiWaypointTrajectory, Scalar> calculate_trajectory_cost_and_optimal_velocities(std::vector<int> &scheduled_locations_idx,
                                                                                             EnvConfig& env_params,
                                                                                             bool sample_start_velocity = false) {

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




constructed_trajectory construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  EnvConfig& env_params,
  MultiWaypointTrajectory mwp_trajectory = MultiWaypointTrajectory{}) {

  // location_positions  &
  // rewards  &
  // t_max  &
  std::vector<Vector<3>>& location_positions = env_params.location_positions;
  std::vector<Scalar>& rewards = env_params.rewards;
  Scalar t_max = env_params.t_max;

  // precalculated_costs  &
  // velocity_norm_samples  &
  // heading_angle_samples  &
  travel_cost_map& precalculated_costs = env_params.precalculated_costs;
  std::vector<Scalar>& velocity_norm_samples = env_params.velocity_norm_samples;
  std::vector<Scalar>& heading_angle_samples = env_params.heading_angle_samples;


  std::vector<Vector<3>> scheduled_locations{};
  for (int idx : scheduled_locations_idx) {
    scheduled_locations.push_back(env_params.location_positions[idx]);
  }
  std::vector<int> unscheduled_locations_idx = get_missing_values_in_range(scheduled_locations_idx,
                                                           scheduled_locations_idx[0],
                                                           scheduled_locations_idx[scheduled_locations_idx.size()-1]);

  MultiWaypointTrajectory current_trajectory;
  if (mwp_trajectory.size() > 0) {
    current_trajectory = mwp_trajectory;
    std::cout << "[NOTE]: trajectory received as param. " << std::endl;
  } else {
    auto current_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
      scheduled_locations_idx,
      env_params,
      true);
    current_trajectory = std::get<0>(current_traj_and_time);
  }


  Scalar current_cost = get_mwp_trajectory_cost(current_trajectory);
  Scalar collected_reward = get_mwp_trajectory_reward(scheduled_locations_idx, rewards);

  while (current_cost < t_max) {

    // <what_idx_to_insert, where_to_insert, velocity>
    std::tuple<int, int, Vector<3>, MultiWaypointTrajectory, Scalar, std::vector<Vector<3>>> best_insertion_so_far{};
    Scalar ratio_of_best_insertion_so_far = -1;
    // Try to schedule every unscheduled location
    for (int unscheduled_idx : unscheduled_locations_idx) {
//      std::cout << unscheduled_idx << std::endl;
      // For every possible insertion_idx
      //                     possible insertion spot
      //                          v
      // Example: {start_location, end_location}
      for (int insertion_idx = 1; insertion_idx < scheduled_locations.size(); insertion_idx++) {
        //        std::cout << "insertion_idx -> " << insertion_idx << std::endl;
        Scalar curr_min_insertion_cost = DBL_MAX;
        int pred_idx = scheduled_locations_idx[insertion_idx-1];
        int succ_idx = scheduled_locations_idx[insertion_idx];

        // For every possible combination of norm and heading_angle
        for (Scalar norm1 : velocity_norm_samples) {
          for (Scalar angle1 : heading_angle_samples) {
            Vector<3> position_to_schedule = location_positions[unscheduled_idx];
            Scalar pred_to_curr_cost, curr_to_succ_cost, pred_to_succ_cost;

            // Predecessor -> Current [cost]
            Scalar pred_norm = current_trajectory[insertion_idx-1].inp_from_v_norm;
            Scalar pred_angle = current_trajectory[insertion_idx-1].inp_from_v_angle;
            pred_to_curr_cost = precalculated_costs[pred_idx][unscheduled_idx][pred_norm][norm1][pred_angle][angle1];
            if (pred_to_curr_cost == 0) {
              QuadState pred_state;
              QuadState curr_state;
              pred_state.setZero();
              curr_state.setZero();
              pred_state.p = location_positions[pred_idx];
              curr_state.p = position_to_schedule;
              pred_state.v = current_trajectory[insertion_idx-1].get_start_state().v;
              curr_state.v = to_velocity_vector(norm1, angle1);
              PointMassTrajectory3D tr(pred_state, curr_state, env_params.max_acc_per_axis, true);
              if (!tr.exists()) {
                std::cout << "Not-existing here -> time -> " << tr.time() << " speed -> " << pred_state.v.transpose() << std::endl;
                continue;
              } else {
                pred_to_curr_cost = tr.time();
              }
            }

            // Current -> Successor [cost]
            Scalar succ_norm = current_trajectory[insertion_idx-1].inp_to_v_norm;
            Scalar succ_angle = current_trajectory[insertion_idx-1].inp_to_v_angle;
            curr_to_succ_cost = precalculated_costs[unscheduled_idx][succ_idx][norm1][succ_norm][angle1][succ_angle];
            if (curr_to_succ_cost == 0) {
              QuadState curr_state;
              QuadState succ_state;
              curr_state.setZero();
              succ_state.setZero();
              curr_state.p = position_to_schedule;
              succ_state.p = location_positions[succ_idx];
              curr_state.v = to_velocity_vector(norm1, angle1);
              succ_state.v = current_trajectory[insertion_idx-1].get_end_state().v;
              PointMassTrajectory3D tr(curr_state, succ_state, env_params.max_acc_per_axis, true);
              if (!tr.exists()) {
                std::cout << "Not-existing here -> time -> " << tr.time() << " speed -> " << succ_state.v.transpose() << std::endl;
                continue;
              } else {
                curr_to_succ_cost = tr.time();
              }
            }

            // Predecessor -> Successor [cost]
            pred_to_succ_cost = current_trajectory[insertion_idx-1].time();
            if (pred_to_curr_cost == 0 || curr_to_succ_cost == 0) {
              std::cout << " This shouldn't happen !!!!!!!!" << std::endl;
            }

            Scalar cost_of_insertion = pred_to_curr_cost + curr_to_succ_cost - pred_to_succ_cost;
            Scalar ratio = rewards[unscheduled_idx] / cost_of_insertion;

            if (ratio > ratio_of_best_insertion_so_far) {
//              std::cout << "true " << std::endl;
              std::vector<Vector<3>> try_scheduled_locations = scheduled_locations;
              std::vector<int> try_scheduled_locations_idx = scheduled_locations_idx;

              try_scheduled_locations.insert(try_scheduled_locations.begin()+insertion_idx, position_to_schedule);
              try_scheduled_locations_idx.insert(try_scheduled_locations_idx.begin()+insertion_idx, unscheduled_idx);

              auto try_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
                try_scheduled_locations_idx,
                env_params,
                true);
              auto new_cost = std::get<1>(try_traj_and_time);
//              Scalar new_cost = cost_of_insertion + current_cost;

              if (new_cost < t_max) {
                // HERE IS WHERE WE LOSE TIME

//                              std::cout << "new cost -> " << new_cost << std::endl;
//                              std::cout << "Real coi -> " << std::get<1>(try_traj_and_time) << std::endl;
                // HERE IS WHERE WE LOSE TIME
                auto new_traj = std::get<0>(try_traj_and_time);
                ratio_of_best_insertion_so_far = ratio;
                Vector<3> vel_vector = to_velocity_vector(norm1, angle1);
                best_insertion_so_far = {unscheduled_idx, insertion_idx, vel_vector, new_traj, new_cost, try_scheduled_locations};
              }
            }
          }
        }
      }
    }
    if (std::get<0>(best_insertion_so_far) == 0 || unscheduled_locations_idx.size() == 0)  {
      // Can't insert any point -> break
      break;
    }

    //    std::cout << "Best ration -> " << ratio_of_best_insertion_so_far << std::endl;
    int what_to = std::get<0>(best_insertion_so_far);
    int where_to = std::get<1>(best_insertion_so_far);

    current_trajectory = std::get<3>(best_insertion_so_far);
    current_cost = std::get<4>(best_insertion_so_far);
    scheduled_locations = std::get<5>(best_insertion_so_far);
    collected_reward += rewards[what_to];

    // Update `scheduled_locations_idx` and Remove scheduled_location from `unscheduled_locations_idx`
    scheduled_locations_idx.insert(scheduled_locations_idx.begin()+where_to, what_to);
    unscheduled_locations_idx.erase(std::remove(unscheduled_locations_idx.begin(), unscheduled_locations_idx.end(), what_to), unscheduled_locations_idx.end());
  }

  std::cout <<  "-------------------------" << std::endl;
  std::cout << "Final cost: " << current_cost << std::endl;
  std::cout << "Collected reward: " << collected_reward << std::endl;
  std::cout << "Actual cost: " << get_mwp_trajectory_cost(current_trajectory) << std::endl;
  std::cout << "Actuual reward: " << get_mwp_trajectory_reward(scheduled_locations_idx, rewards) << std::endl;

  // OUTPUT: MultiWaypointTrajectory, total_cost, total_reward, scheduled_locations_idx, unscheduled_locations_idx
  return {current_trajectory, current_cost, collected_reward, scheduled_locations_idx, unscheduled_locations_idx};
}
void destruction_heuristic_3(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             std::vector<Scalar> ratios,
                             EnvConfig& env_state_config) {

  int min_idx = 1;
  Scalar min_heu = 10000000000;
  for (int i = 1; i < sched_loc.size()-1; i++) {

    std::vector<int> to_try{sched_loc[i-1], sched_loc[i], sched_loc[i+1]};
    auto to_try_optimal = calculate_trajectory_cost_and_optimal_velocities(to_try,
                                                                           env_state_config,
                                                                           true);
    Scalar optimal_time = std::get<1>(to_try_optimal);
    Scalar curr_time = curr_traj[i-1].time() + curr_traj[i].time();
    Scalar diff = abs(optimal_time - curr_time);
    Scalar heu = ratios[i] / diff;
    if (heu < min_heu) {
      min_heu = heu;
      min_idx = i;
    }
  }
  unsched_loc.push_back(sched_loc[min_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[min_idx]), sched_loc.end());
}

void destruction_heuristic_2(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             EnvConfig& env_state_config) {

  int max_idx = 1;
  Scalar max_diff = -1;

  for (int i = 1; i < sched_loc.size()-1; i++) {
    std::vector<int> to_try{sched_loc[i-1], sched_loc[i], sched_loc[i+1]};
    auto to_try_optimal = calculate_trajectory_cost_and_optimal_velocities(to_try,
                                                                           env_state_config,
                                                                           true);
    Scalar optimal_time = std::get<1>(to_try_optimal);
    Scalar curr_time = curr_traj[i-1].time() + curr_traj[i].time();
    Scalar diff = abs(optimal_time - curr_time);
    if (diff < max_diff) {
      max_diff = diff;
      max_idx = i;
    }
//    std::cout << "Time now -> " << mvt[i-1].time() + mvt[i].time() << " | time Optimal -> " << std::get<1>(to_try_optimal) << std::endl;
//    std::cout << "Diff -> " << abs(mvt[i-1].time() + mvt[i].time() - std::get<1>(to_try_optimal)) << std::endl;

  }
  unsched_loc.push_back(sched_loc[max_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[max_idx]), sched_loc.end());
}


void destruction_heuristic_1(std::vector<int>& sched_loc,
                            std::vector<int>& unsched_loc,
                            std::vector<Scalar>& ratios) {
//  for (auto asas : unsched_loc) {
//    std::cout  << asas << " -> ";
//  }
//  std::cout << std::endl;
  int min_idx = 1;
  Scalar curr_min_ratio = ratios[min_idx];
  for (int i = 1; i < sched_loc.size()-1; i++) {
    if (ratios[i] < curr_min_ratio) {
      curr_min_ratio = ratios[i];
      min_idx = i;
    }
  }
//  std::cout << "To remove loc -> " << sched_loc[min_idx] << " with ration of " << ratios[min_idx] << std::endl;
  unsched_loc.push_back(sched_loc[min_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[min_idx]), sched_loc.end());
//  std::cout<< "Push " << sched_loc[min_idx] << std::endl;


}

std::vector<Scalar> calculate_heuristic_ratio(std::vector<int>& scheduled_locations_idx,
                                              MultiWaypointTrajectory& current_trajectory,
                                              std::vector<Scalar> rewards,
                                              travel_cost_map& precalculated_costs) {

  std::vector<Scalar> final_ratios{-1};
  // Finding final ratios
  for (int i = 1; i < scheduled_locations_idx.size()-1; i++) {
    Scalar pred_to_curr_cost, curr_to_succ_cost, pred_to_succ_cost;
    Scalar pred_norm = current_trajectory[i-1].inp_from_v_norm;
    Scalar pred_angle = current_trajectory[i-1].inp_from_v_angle;
    Scalar curr_norm = current_trajectory[i-1].inp_to_v_norm;
    Scalar curr_angle = current_trajectory[i-1].inp_to_v_angle;
    Scalar succ_norm = current_trajectory[i].inp_from_v_norm;
    Scalar succ_angle = current_trajectory[i].inp_from_v_angle;

    // Predecessor -> Current [cost]
    pred_to_curr_cost = precalculated_costs[scheduled_locations_idx[i-1]][scheduled_locations_idx[i]][pred_norm][curr_norm][pred_angle][curr_angle];
    // Current -> Successor [cost]
    curr_to_succ_cost = precalculated_costs[scheduled_locations_idx[i]][scheduled_locations_idx[i+1]][curr_norm][succ_norm][curr_angle][succ_angle];
    // Predecessor -> Successor [cost]
    pred_to_succ_cost = precalculated_costs[scheduled_locations_idx[i-1]][scheduled_locations_idx[i+1]][pred_norm][succ_norm][pred_angle][succ_angle];

    Scalar cost_of_insertion = pred_to_curr_cost + curr_to_succ_cost - pred_to_succ_cost;
    Scalar ratio = rewards[scheduled_locations_idx[i]] / cost_of_insertion;
    final_ratios.push_back(ratio);
    //    std::cout << " ratio -> " << ratio << std::endl;
    //    if (ratio < 5) {
    //      std::cout << "reward -> " << rewards[scheduled_locations_idx[i]] <<" " << scheduled_locations_idx[i] <<  std::endl;
    //    }
  }
  final_ratios.push_back(-1);
  //  std::cout  << std::endl;
  //  for (auto tr : current_trajectory) {
  //    std::cout << "from -> " << tr.get_start_state().p.transpose() << " to -> " << tr.get_end_state().p.transpose() << std::endl;
  //  }
  return final_ratios;
}


std::vector<int> destruction_heuristic_paper(const constructed_trajectory& constr_tr,
                                 Scalar percentage,
                                 EnvConfig& env_params) {
  // ENV PARAMS
  travel_cost_map& travel_costs = env_params.precalculated_costs;
  std::vector<Scalar>& rewards = env_params.rewards;

  MultiWaypointTrajectory mvt = std::get<0>(constr_tr);
  Scalar cost = std::get<1>(constr_tr);
  Scalar reward = std::get<2>(constr_tr);
  std::vector<int> sched_loc = std::get<3>(constr_tr);
  std::vector<int> unsched_loc = std::get<4>(constr_tr);
  std::vector<Scalar> ratios = calculate_heuristic_ratio(sched_loc, mvt, rewards, travel_costs);

  int positions_to_remove = sched_loc.size() * percentage / 100;

//  std::cout << "Percentage -> " << percentage <<"% " << "removes " << positions_to_remove << " positions" << std::endl;

  for (int i = 0; i < positions_to_remove; i++) {


    int random_number = randint(1, 3);
//    std::cout << random_number << std::endl;
    if (random_number == 1) {
      destruction_heuristic_1(sched_loc, unsched_loc, ratios);
    } else if (random_number == 2) {
      destruction_heuristic_2(sched_loc, unsched_loc, mvt, env_params);
    } else if (random_number == 3) {
      destruction_heuristic_3(sched_loc, unsched_loc, mvt, ratios, env_params);
    }
    auto new_traj_time = calculate_trajectory_cost_and_optimal_velocities(sched_loc,env_params, true);
    mvt = std::get<0>(new_traj_time);
    cost = std::get<1>(new_traj_time);
//    std::cout << cost << std::endl;
    ratios = calculate_heuristic_ratio(sched_loc, mvt, rewards, travel_costs);
  }
//  std::cout << "Final cost " << cost << std::endl;
//  for (auto sl : sched_loc) {
//    std::cout  << sl << " -> ";
//  }
//  std::cout << std::endl;
//  for (auto asas : unsched_loc) {
//    std::cout  << asas << " -> ";
//  }
  return sched_loc;
}


constructed_trajectory run_paper_heuristic(EnvConfig& env_state_config) {

  // idx == location by idx in `location_positions`
  std::vector<int> scheduled_locations_idx = {0, (int)env_state_config.location_positions.size()-1};

  constructed_trajectory initial_constr = construction_heuristic(
    scheduled_locations_idx,
    env_state_config);

  constructed_trajectory best_constr_yet = initial_constr;
  MultiWaypointTrajectory best_tr_yet = std::get<0>(initial_constr);
  Scalar best_cost = std::get<1>(initial_constr);
  Scalar best_reward_yet = std::get<2>(initial_constr);
  std::vector<int> best_scheduled_positions = std::get<3>(initial_constr);


  for (int j = 0; j < 100; j++) {
    std::cout << "First Construction #" << j << std::endl;

    scheduled_locations_idx = destruction_heuristic_paper(initial_constr, 50, env_state_config);
    initial_constr = construction_heuristic(scheduled_locations_idx, env_state_config);
    if (std::get<2>(initial_constr) > best_reward_yet) {
      best_constr_yet = initial_constr;
      best_tr_yet = std::get<0>(initial_constr);
      best_cost = std::get<1>(initial_constr);
      best_reward_yet = std::get<2>(initial_constr);
      best_scheduled_positions = std::get<3>(initial_constr);
    }
  }

//  for (int j = 0; j < 100; j++) {
//    std::cout << "Second Construction #" << j << std::endl;
//    auto destroyed_solution = destruction_heuristic_paper(initial_constr, 20, env_state_config);
//    scheduled_locations_idx = std::get<0>(destroyed_solution);
//    initial_constr = construction_heuristic(scheduled_locations_idx, env_state_config);
//    if (std::get<2>(initial_constr) > best_reward_yet) {
//      best_constr_yet = initial_constr;
//      best_tr_yet = std::get<0>(initial_constr);
//      best_cost = std::get<1>(initial_constr);
//      best_reward_yet = std::get<2>(initial_constr);
//      best_scheduled_positions = std::get<3>(initial_constr);
//    }
//  }

//  std::cout << "------------------  STATS  ---------------" << std::endl;
//  std::cout << "Final reward -> " << best_reward_yet << std::endl;
//  std::cout << "Final cost -> " << best_cost << std::endl;
//  std::cout << std::endl;
//  std::cout << "Velocity norms:" << std::endl;
//  auto vel_norms = get_mwp_trajectory_velocities(best_tr_yet);
//  for (auto n : vel_norms) {
//    std::cout << n << ", ";
//  }
//  std::cout << std::endl;
//  std::cout << "Velocity angles:" << std::endl;
//  auto yaw_angles = get_mwp_trajectory_yaw_angles(best_tr_yet);
//  for (auto ya : yaw_angles) {
//    std::cout << ya << ", ";
//  }
//  std::cout << std::endl;
//  std::cout << "Scheduled locations:" << std::endl;
//  for (auto sch: best_scheduled_positions) {
//    std::cout << sch << " -> ";
//  }
//  std::cout << std::endl;

  return best_constr_yet;
}

void get_positions_travel_costs(std::string config_file, int argc, char** cli_args)
{
  // LOAD CONFIG
  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();

  exit(1);

  auto heu_result = run_paper_heuristic(env_state_config);
  MultiWaypointTrajectory best_tr_yet = std::get<0>(heu_result);
  Scalar best_cost = std::get<1>(heu_result);
  Scalar best_reward_yet = std::get<2>(heu_result);
  std::vector<int> best_scheduled_positions = std::get<3>(heu_result);

  std::cout << "------------------  FINAL STATS AFTER RUN PAPER HEU  ---------------" << std::endl;
  std::cout << "Reward: " << best_reward_yet << std::endl;
  std::cout << "Cost: " << best_cost << std::endl;
  std::cout << "Scheduled locations:" << std::endl;
  for (auto sch: best_scheduled_positions) {
    std::cout << sch << " -> ";
  }
  VelocitySearchGraph::saveTrajectoryEquitemporal(best_tr_yet, "samples_pmm.csv");
  std::cout <<  "--------------------" << std::endl;
  std::cout << "Saved equitemporal." << std::endl;
  exit(1);

  std::cout << "------------------  START CONE REFOCUSING  ---------------" << std::endl;
  std::vector<Scalar> vel_norms = get_mwp_trajectory_velocities(best_tr_yet);
  std::vector<Scalar> yaw_angles = get_mwp_trajectory_yaw_angles(best_tr_yet);
  std::vector<Vector<3>> scheduled_loc_gates;
  for (int pos_idx : best_scheduled_positions) {
    scheduled_loc_gates.push_back(env_state_config.location_positions[pos_idx]);
  }
  Vector<3> start_vel = to_velocity_vector(best_tr_yet[0].inp_from_v_norm, best_tr_yet[0].inp_from_v_angle);

  auto cone_refoces_tr = test_pmm(scheduled_loc_gates, vel_norms, yaw_angles, start_vel);


  std::cout << "------------------  AUGMENTING  ---------------" << std::endl;
  constructed_trajectory augmented_cone_refocused_tr = construction_heuristic(
    best_scheduled_positions,
    env_state_config,
    cone_refoces_tr);

  auto tr = std::get<0>(augmented_cone_refocused_tr);
  auto cost = std::get<1>(augmented_cone_refocused_tr);
  auto reward = std::get<2>(augmented_cone_refocused_tr);
  auto scheduled_idx = std::get<3>(augmented_cone_refocused_tr);
  Scalar actual_cost = get_mwp_trajectory_cost(tr);
  Scalar actual_reward = get_mwp_trajectory_reward(scheduled_idx, env_state_config.rewards);
  std::cout << " cost -> " << cost << std::endl;
  std::cout << " actual cost -> " << actual_cost << std::endl;
  std::cout << " reward -> " << reward << std::endl;
  std::cout << " actual reward -> " << actual_reward << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(tr, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;

  std::vector<Scalar> vel_norms2 = get_mwp_trajectory_velocities(tr);
  std::vector<Scalar> yaw_angles2 = get_mwp_trajectory_yaw_angles(tr);
  std::vector<Vector<3>> scheduled_loc_gates2;
  for (int pos_idx : scheduled_idx) {
    scheduled_loc_gates2.push_back(env_state_config.location_positions[pos_idx]);
  }
  Vector<3> start_vel2 = to_velocity_vector(tr[0].inp_from_v_norm, tr[0].inp_from_v_angle);
  auto cone_refoces_tr1 = test_pmm(scheduled_loc_gates2, vel_norms2, yaw_angles2, start_vel2);
  Scalar final_actual_cost = get_mwp_trajectory_cost(cone_refoces_tr1);
  Scalar final_actual_reward = get_mwp_trajectory_reward(scheduled_idx, env_state_config.rewards);
  std::cout << " Final  cost -> " << final_actual_cost << std::endl;
  std::cout << " Final reward -> " << final_actual_reward << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(cone_refoces_tr1, "samples_pmm.csv");
  std::cout << "Saved final equitemporal." << std::endl;

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2 22222222222222

  constructed_trajectory augmented_cone_refocused_tr2 = construction_heuristic(
    scheduled_idx,
    env_state_config,
    cone_refoces_tr1);

  auto tr3 = std::get<0>(augmented_cone_refocused_tr2);
  auto cost3 = std::get<1>(augmented_cone_refocused_tr2);
  auto reward3 = std::get<2>(augmented_cone_refocused_tr2);
  auto scheduled_idx3 = std::get<3>(augmented_cone_refocused_tr2);
  Scalar actual_cost3 = get_mwp_trajectory_cost(tr3);
  Scalar actual_reward3 = get_mwp_trajectory_reward(scheduled_idx3, env_state_config.rewards);
  std::cout << " cost -> " << cost3 << std::endl;
  std::cout << " actual cost -> " << actual_cost3 << std::endl;
  std::cout << " reward -> " << reward3 << std::endl;
  std::cout << " actual reward -> " << actual_reward3 << std::endl;


  std::vector<Scalar> vel_norms3 = get_mwp_trajectory_velocities(tr3);
  std::vector<Scalar> yaw_angles3 = get_mwp_trajectory_yaw_angles(tr3);
  std::vector<Vector<3>> scheduled_loc_gates3;
  for (int pos_idx : scheduled_idx3) {
    scheduled_loc_gates3.push_back(env_state_config.location_positions[pos_idx]);
  }
  Vector<3> start_vel3 = to_velocity_vector(tr3[0].inp_from_v_norm, tr3[0].inp_from_v_angle);

  auto cone_refoces_tr3 = test_pmm(scheduled_loc_gates3, vel_norms3, yaw_angles3, start_vel3);
  VelocitySearchGraph::saveTrajectoryEquitemporal(tr3, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;
  std::cout << "Cost -> " << get_mwp_trajectory_cost(cone_refoces_tr3) << std::endl;
  std::cout << "Reward -> " << get_mwp_trajectory_reward(scheduled_idx3, env_state_config.rewards) << std::endl;
}

int main(int argc, char** argv) {
//  std::cout << to_velocity_vector(1.06066, 1.58) << std::endl;
//  exit(1);
//  test_pmm();
//
  srand(13);
  get_positions_travel_costs("/Users/markv/pmm_planner/new_config.yaml", argc, argv);

  return 0;
}
