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


// point mass model
// Keep this code for code refocusing example
int test_pmm(int argc, char** argv) {
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
  // Maximal values of the samples.
  const Scalar yaw_pitch_cone_angle_boundary = 60.0;

  // Velocity norm samples.
  const Scalar min_velocity_norm = 5.0;
  const Scalar min_velocity_norm_boundary = 1.0;
  const Scalar max_velocity_norm = 17.0;
  const Scalar precision_velocity_norm = 8.0;
  const Scalar max_acc_norm = 32.94;
  const bool end_free = config["end_free"].as<bool>();
  Scalar max_yaw_pitch_ang = config["max_yaw_pitch_ang"].as<Scalar>();
  Scalar precision_yaw_pitch_ang = config["precision_yaw_pitch_ang"].as<Scalar>();

  // Construct a graph given all the parameters of the environment.
  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  Vector<3> start_velocity;
  Vector<3> start_position;
  Vector<3> end_velocity;
  Vector<3> end_position;
  config["start"]["velocity"] >> start_velocity;
  config["start"]["position"] >> start_position;
  config["end"]["velocity"] >> end_velocity;
  config["end"]["position"] >> end_position;

  // Vector of 'heading angles'(positions) of all gates (including start and end).
  std::vector<Vector<3>> gates_waypoints;
  if (!parseArrayParam<Vector<3>>(config, "locations", gates_waypoints))
    std::cerr << "can't load param gates" << std::endl;
  gates_waypoints.insert(gates_waypoints.begin(), start_position);
  gates_waypoints.push_back(end_position);
  std::vector<Scalar> gates_yaw_deg;
  if (!parseArrayParam<Scalar>(config, "gates_orientations", gates_yaw_deg))
    std::cerr << "can't load param gates_orientations" << std::endl;

  std::vector<Scalar> gates_pitch_deg;
  if (!parseArrayParam<Scalar>(config, "gates_pitch_orientations",
                               gates_pitch_deg))
    std::cerr << "can't load param gates_pitch_orientations" << std::endl;

  std::vector<Scalar> gates_vel_norms;
  gates_vel_norms.resize(gates_pitch_deg.size(),
                         (max_velocity_norm + min_velocity_norm) / 2.0);



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

  for (size_t i = 0; i < tr.size(); i++) {
    const Vector<3> vel = tr[i].get_end_state().v;
    const Vector<3> vel_norm = vel.normalized();
    const Scalar pitch = asin(-vel_norm(2)) * 180 / M_PI;
    const Scalar yaw = atan2(vel_norm(1), vel_norm(0)) * 180 / M_PI;
    gates_yaw_deg[i + 1] = yaw;
    gates_pitch_deg[i + 1] = pitch;
    gates_vel_norms[i + 1] = vel.norm();
  }
  find_vel.toc();
  find_vel.print();

//  std::cout << "output tr size " << tr.size() << std::endl;
  for (size_t i = 0; i < tr.size(); i++) {
    if (i == 0) {
      Vector<3> vel_norm = tr[i].get_start_state().v.normalized();
      Scalar pitch = asin(-vel_norm(2)) * 180 / M_PI;
      Scalar yaw = atan2(vel_norm(1), vel_norm(0)) * 180 / M_PI;
//      std::cout << i << " pos " << tr[i].get_start_state().p.transpose()
//                << " vel " << tr[i].get_start_state().v.transpose() << " acc "
//                << tr[i].get_start_state().a.transpose() << " thrust norm "
//                << (tr[i].get_start_state().a - GVEC).norm() << " exists "
//                << tr[i].exists() << " yaw " << yaw << " pitch " << pitch
//                << std::endl;
    }
    Vector<3> vel_norm = tr[i].get_end_state().v.normalized();
    Scalar pitch = asin(-vel_norm(2)) * 180 / M_PI;
    Scalar yaw = atan2(vel_norm(1), vel_norm(0)) * 180 / M_PI;
//    std::cout << i + 1 << " pos " << tr[i].get_end_state().p.transpose()
//              << " vel " << tr[i].get_end_state().v.transpose() << " acc "
//              << tr[i].get_end_state().a.transpose() << " thrust norm "
//              << (tr[i].get_end_state().a - GVEC).norm() << " exists "
//              << tr[i].exists() << " yaw " << yaw << " pitch " << pitch
//              << std::endl;

    sum_times += tr[i].time();
    if (tr[i].time() - tr[i].time_min() > PRECISION_PMM_VALUES) {
//      std::cout << "bad time!!!!!!" << std::endl;
    }
  }
  std::cout <<  "--------------------" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(tr, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(tr,
                                                 "samples_equidistant.csv");
  std::cout << "Saved equidistant." << std::endl;

  return 0;
}



std::tuple<MultiWaypointTrajectory, Scalar> calculate_trajectory_cost_and_optimal_velocities(
                                                                                             std::vector<int> &scheduled_locations_idx,
                                                                                             EnvConfig& env_params) {
  travel_cost_map& travel_costs = env_params.precalculated_costs;
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>>& velocity_sample_tuples = env_params.velocity_samples_tuples;
  Vector<3> start_velocity = env_params.start_velocity;
//  Vector<3> &end_velocity = env_params.location_positions[scheduled_locations_idx[scheduled_locations_idx.size()-1]];
  Vector<3> &max_acc_per_axis = env_params.max_acc_per_axis;
  std::vector<Vector<3>> location_positions{};

  for (int idx : scheduled_locations_idx) {
    location_positions.push_back(env_params.location_positions[idx]);
//    if (scheduled_locations_idx.size() == 3) {
//      std::cout<< idx << std::endl;
//    }
  }
//  if (scheduled_locations_idx.size() == 3) {
//    std::cout << std::endl;
//  }



  Scalar shortest_time = DBL_MAX;
  int number_of_locations = location_positions.size();
  int number_of_velocity_samples = velocity_sample_tuples.size();
  std::tuple<Scalar, Scalar> start_vel_mag_ang = to_velocity_norm_and_angle(start_velocity);

  // {id_of_sample_in_location_before, time_from_the_start}, idx - location_idx
  std::vector<std::vector<std::pair<int, Scalar>>> shortest_samples_times;
  shortest_samples_times.resize(number_of_locations);

  // For 'start' location, it takes 0 seconds to reach.
  shortest_samples_times[0].push_back({-1, 0});
  const std::pair<int, Scalar> not_reached = std::make_pair(-1, DBL_MAX);
  for (int loc_id = 1; loc_id < number_of_locations; loc_id++) {
    shortest_samples_times[loc_id].resize(number_of_velocity_samples, {-1, DBL_MAX});
    std::fill(shortest_samples_times[loc_id].begin(),
              shortest_samples_times[loc_id].end(), not_reached);
  }

  // Samples have [loc_id][sample_id][{velocity, norm, heading_angle)}]
  std::vector<std::vector<std::tuple<Vector<3>, Scalar, Scalar>>> gate_velocity_samples;
  gate_velocity_samples.resize(number_of_locations);
  for (int loc_id = 1; loc_id < number_of_locations; loc_id++) {
    gate_velocity_samples[loc_id] = velocity_sample_tuples;
  }
  gate_velocity_samples[0].push_back({start_velocity, 0, 0});


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
    for (int idx_from = 0; idx_from < from_shortest_samples_dists.size(); idx_from++) {
      const Scalar& time_from = from_shortest_samples_dists[idx_from].second;
      const Vector<3>& from_v = std::get<0>(from_samples[idx_from]);
      if (time_from == DBL_MAX) {
        // No path reaching that sample.
        continue;
      }
      // loop for indexes of the samples in the gate gid_to
      for (int idx_to = 0; idx_to < to_shortest_samples_dists.size(); idx_to++) {
        std::pair<int, Scalar>& shortest_to = to_shortest_samples_dists[idx_to];

        Scalar from_norm = std::get<1>(from_samples[idx_from]);
        Scalar from_angle = std::get<2>(from_samples[idx_from]);
        Scalar to_norm = std::get<1>(to_samples[idx_to]);
        Scalar to_angle = std::get<2>(to_samples[idx_to]);

        // If location is START, the velocity is fixed (not sampled).
        if (loc_to-1 == 0) {
          from_norm = std::get<0>(start_vel_mag_ang);
          from_angle = std::get<1>(start_vel_mag_ang);
        }
        if (loc_to == 0) {
          to_norm = std::get<0>(start_vel_mag_ang);
          to_angle = std::get<1>(start_vel_mag_ang);
        }

        Scalar time_between = travel_costs[scheduled_locations_idx[loc_to-1]][scheduled_locations_idx[loc_to]][from_norm][to_norm][from_angle][to_angle];
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
  //  std::cout << "shortest_time " << shortest_time << std::endl << std::endl;

  // Knowing the id of the sample of the last location that gives us the best
  // time, we can traceback the optimal velocities at all previous locations.


  std::vector<std::tuple<Vector<3>, Scalar, Scalar>> found_optimal_gates_velocities;
  // We don't find optimal velocities for 'start', thus, number_of_locations - 1.
  found_optimal_gates_velocities.resize(number_of_locations - 1);
  int prev_sample_idx = end_best_idx;
  for (int loc_id = endi; loc_id > 0; loc_id--) {
    found_optimal_gates_velocities[loc_id - 1] = gate_velocity_samples[loc_id][prev_sample_idx];
    prev_sample_idx = shortest_samples_times[loc_id][prev_sample_idx].first;
  }

  // Populate trajectories between all points.
  MultiWaypointTrajectory trajectories;
  // Total time of all trajectories.
  Scalar time_sum = 0;

  // ------ Add trajectory "start" -> "location_1"
  QuadState from_start_state;
  from_start_state.p = location_positions[0];
  from_start_state.v = start_velocity;
  from_start_state.v_norm = std::get<0>(start_vel_mag_ang);
  from_start_state.v_angle = std::get<1>(start_vel_mag_ang);;
  QuadState to_state;
  to_state.p = location_positions[1];
  to_state.v = std::get<0>(found_optimal_gates_velocities[0]);
  to_state.v_norm = std::get<1>(found_optimal_gates_velocities[0]);
  to_state.v_angle = std::get<2>(found_optimal_gates_velocities[0]);
  PointMassTrajectory3D tr_max_acc_start;
  tr_max_acc_start =
    PointMassTrajectory3D(from_start_state, to_state, max_acc_per_axis, true);

  trajectories.push_back(tr_max_acc_start);
  time_sum += tr_max_acc_start.time();
  if (!tr_max_acc_start.exists()) {
    std::cout << "Trajectory does not exist from start location." << std::endl;
    return {MultiWaypointTrajectory(), DBL_MAX};
  }
  // ------ Add trajectory "start" -> "location_1" ---- END

  // ------ Add all other trajectories
  for (size_t i = 2; i < location_positions.size(); i++) {
    QuadState from_state_b;
    from_state_b.p = location_positions[i - 1];
    from_state_b.v = std::get<0>(found_optimal_gates_velocities[i - 2]);
    from_state_b.v_norm = std::get<1>(found_optimal_gates_velocities[i - 2]);
    from_state_b.v_angle = std::get<2>(found_optimal_gates_velocities[i - 2]);

    QuadState to_state_b;
    to_state_b.p = location_positions[i];
    to_state_b.v = std::get<0>(found_optimal_gates_velocities[i - 1]);
    to_state_b.v_norm = std::get<1>(found_optimal_gates_velocities[i - 1]);
    to_state_b.v_angle = std::get<2>(found_optimal_gates_velocities[i - 1]);

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
  //  std::cout << "End of trajectory calculation. Return trajectories found.." << std::endl;
  return {trajectories, shortest_time};
}




std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>, std::vector<int>> construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  std::vector<int> unscheduled_locations_idx,
  EnvConfig& env_params) {
  // CFG
  // location_positions  &
  // start_velocity  &
  // end_velocity  &
  // max_acc_per_axis  &
  // t_max  &
  // rewards  &
  std::vector<Vector<3>>& location_positions = env_params.location_positions;
  Vector<3>& start_velocity = env_params.start_velocity;
  Vector<3>& end_velocity = env_params.end_velocity;
  Vector<3>& max_acc_per_axis = env_params.max_acc_per_axis;
  Scalar t_max = env_params.t_max;
  std::vector<Scalar>& rewards = env_params.rewards;


  // PRECALC
  // precalculated_costs  &
  // velocity_samples_tuples  &
  // velocity_norm_samples  &
  // heading_angle_samples  &
  travel_cost_map& precalculated_costs = env_params.precalculated_costs;
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>>& velocity_samples_tuples = env_params.velocity_samples_tuples;
  std::vector<Scalar>& velocity_norm_samples = env_params.velocity_norm_samples;
  std::vector<Scalar>& heading_angle_samples = env_params.heading_angle_samples;

  // unscheduled_locations_idx
  // scheduled_locations_idx
  // scheduled_locations -> get from `location_positions` by `scheduled_locations_idx`
  std::vector<Vector<3>> scheduled_locations{};
  for (int idx : scheduled_locations_idx) {
    scheduled_locations.push_back(env_params.location_positions[idx]);
  }


  auto current_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
    scheduled_locations_idx,
    env_params);
  MultiWaypointTrajectory current_trajectory = std::get<0>(current_traj_and_time);
  Scalar current_cost = std::get<1>(current_traj_and_time);
  Scalar collected_reward = 0;
  std::cout << "Current cost " << current_cost << std::endl;

  while (current_cost < t_max) {
    // <what_idx_to_insert, where_to_insert, velocity>
    std::tuple<int, int, Vector<3>, MultiWaypointTrajectory, Scalar, std::vector<Vector<3>>> best_insertion_so_far{};
    Scalar ratio_of_best_insertion_so_far = -1;
    // Try to schedule every unscheduled location
    for (int unscheduled_idx : unscheduled_locations_idx) {
      //      std::cout << "unschedule index -> " << unscheduled_idx << std::endl;
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

            // Current -> Successor [cost]
            Scalar succ_norm = current_trajectory[insertion_idx-1].inp_to_v_norm;
            Scalar succ_angle = current_trajectory[insertion_idx-1].inp_to_v_angle;
            curr_to_succ_cost = precalculated_costs[unscheduled_idx][succ_idx][norm1][succ_norm][angle1][succ_angle];

            // Predecessor -> Successor [cost]
            pred_to_succ_cost = current_trajectory[insertion_idx-1].time();

            Scalar cost_of_insertion = pred_to_curr_cost + curr_to_succ_cost - pred_to_succ_cost;
            Scalar ratio = rewards[unscheduled_idx] / cost_of_insertion;

            if (ratio > ratio_of_best_insertion_so_far) {
              std::vector<Vector<3>> try_scheduled_locations = scheduled_locations;
              std::vector<int> try_scheduled_locations_idx = scheduled_locations_idx;

              try_scheduled_locations.insert(try_scheduled_locations.begin()+insertion_idx, position_to_schedule);
              try_scheduled_locations_idx.insert(try_scheduled_locations_idx.begin()+insertion_idx, unscheduled_idx);

              auto try_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
                                                                                        try_scheduled_locations_idx,
                                                                                        env_params);
              auto new_traj = std::get<0>(try_traj_and_time);
              auto new_cost = std::get<1>(try_traj_and_time);

              if (new_cost < t_max) {
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
    Vector<3> velocity = std::get<2>(best_insertion_so_far);
    MultiWaypointTrajectory new_trajectory = std::get<3>(best_insertion_so_far);
    Scalar new_cost = std::get<4>(best_insertion_so_far);
    std::vector<Vector<3>> new_scheduled_locations = std::get<5>(best_insertion_so_far);
    //    std::cout << "Best ration act -> " << "Insert " << what_to << " at " << where_to << " with velocity ->" << velocity.transpose() <<  std::endl;
    //    std::cout << "New cost -> " << new_cost << " at " << where_to << " with velocity ->" << velocity.transpose() <<  std::endl;

    scheduled_locations = new_scheduled_locations;
    //    std::cout << "size_new_sched_loc " << scheduled_locations.size() << std::endl;
    current_trajectory = new_trajectory;
    current_cost = new_cost;
    collected_reward += rewards[what_to];
    scheduled_locations_idx.insert(scheduled_locations_idx.begin()+where_to, what_to);
    unscheduled_locations_idx.erase(std::remove(unscheduled_locations_idx.begin(), unscheduled_locations_idx.end(), what_to), unscheduled_locations_idx.end());
  }

  std::cout <<  "-------------------------" << std::endl;
  std::cout << "Current cost: " << current_cost << std::endl;
  std::cout << "Collected reward: " << collected_reward << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(current_trajectory, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;
//  for (auto ul : unscheduled_locations_idx) {
//    std::cout << "Unscheduled location -> " << ul << std::endl;
//  }
  for (auto sl : scheduled_locations_idx) {
    std::cout  << sl << " -> ";
  }
//  std::cout  << std::endl;
//  for (auto tr : current_trajectory) {
//    std::cout << "from -> " << tr.get_start_state().p.transpose() << " to -> " << tr.get_end_state().p.transpose() << std::endl;
//  }
  // OUTPUT: MultiWaypointTrajectory, total_cost, total_reward, scheduled_locations_idx, unscheduled_locations_idx
  return {current_trajectory, current_cost, collected_reward, scheduled_locations_idx, unscheduled_locations_idx};
}


void get_positions_travel_costs(std::string config_file)
{
  // LOAD CONFIG
//  std::string cfg_file = config_file;
//  YAML::Node config = YAML::LoadFile(cfg_file);

  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();

//  Vector<3> start_velocity;
//  Vector<3> start_position;
//  Vector<3> end_velocity;
//  Vector<3> end_position;
//  config["start"]["velocity"] >> start_velocity;
//  config["start"]["position"] >> start_position;
//  config["end"]["velocity"] >> end_velocity;
//  config["end"]["position"] >> end_position;
//
//  int V = config["number_of_velocity_samples"].as<int>();
//  int H = config["number_of_angle_samples"].as<int>();
//  Scalar t_max = config["t_max"].as<Scalar>();
//  Scalar max_velocity = config["max_velocity_norm"].as<Scalar>();
//  Scalar max_acceleration = config["max_acceleration"].as<Scalar>();
//  Vector<3> max_acc_per_axis = Vector<3>::Constant(max_acceleration);
//
//  // Vector of 'heading angles'(positions) of all gates (including start and end).
//  std::vector<Vector<3>> location_positions;
//  if (!parseArrayParam<Vector<3>>(config, "locations", location_positions)) {
//    std::cerr << "YAML:PARSER -> Can't load param 'locations'" << std::endl;
//  }
//  location_positions.insert(location_positions.begin(), start_position);
//  location_positions.push_back(end_position);
//
//  std::vector<Scalar> rewards{};
//  if (!parseArrayParam<Scalar>(config, "rewards", rewards)) {
//    std::cerr << "Can't load param 'rewards'" << std::endl;
//  }
//  rewards.insert(rewards.begin(), 0);
//  rewards.push_back(0);
//  // LOAD CONFIG END -------------------------------------------------------
//
//  std::vector<Scalar> velocity_norm_samples = get_velocity_norm_samples(V, max_velocity);
//  std::vector<Scalar> heading_angle_samples = get_heading_angle_samples(H);
//
//  // Samples combinations in vector form. [norm][angle] -> *vector [x, y, z]*
//  norm_angle_to_velocity_vector_map norm_angle_to_vector{};
//  populate_norm_angle_to_velocity_vector_map(norm_angle_to_vector,
//                                             velocity_norm_samples,
//                                             heading_angle_samples);
//
//  // Precalculated costs map.
//  travel_cost_map precalculated_costs{};
//  populate_precalculated_travel_costs_map(precalculated_costs,
//                                          location_positions,
//                                          velocity_norm_samples,
//                                          heading_angle_samples,
//                                          norm_angle_to_vector,
//                                          max_acc_per_axis,
//                                          start_velocity,
//                                          end_velocity);
//
//  // Tuples of all samples <velocity_vector, norm, angle>
//  std::vector<std::tuple<Vector<3>, Scalar, Scalar>> velocity_samples_tuples =
//    generate_total_sample_tuples(norm_angle_to_vector,
//                                 velocity_norm_samples,
//                                 heading_angle_samples);
//  std::vector<int> location_idx{};
//  for (int idx = 0; idx < location_positions.size(); idx++) {
//    location_idx.push_back(idx);
//  }
//  auto traj_and_time_1 = calculate_trajectory_cost_and_optimal_velocities(precalculated_costs,
//                                                   location_positions,
//                                                   velocity_samples_tuples,
//                                                   start_velocity,
//                                                   end_velocity,
//                                                   max_acc_per_axis,
//                                                   location_idx);
//
//  MultiWaypointTrajectory traj = std::get<0>(traj_and_time_1);
//  Scalar traj_time = std::get<1>(traj_and_time_1);
//  std::cout << "Final time is " << traj_time << std::endl;
//  std::cout << "Size " << traj.size() << std::endl;
//  exit(1);
////  std::cout << "Size " << traj[0].inp_from_v_norm << std::endl;
//  VelocitySearchGraph::saveTrajectoryEquitemporal(traj, "samples_pmm.csv");
//  std::cout << "Saved equitemporal." << std::endl;
//  VelocitySearchGraph::saveTrajectoryEquidistant(traj, "samples_equidistant.csv");
//  exit(1);
  // ----------------------------- HEURISTIC IMPLEMENTATION

  // idx == location by idx in `location_positions`
  std::vector<int> scheduled_locations_idx = {0, (int)env_state_config.location_positions.size()-1};
  std::vector<int> unscheduled_locations_idx{};
  for (int i = 1; i < env_state_config.location_positions.size() - 1; i++) {
    unscheduled_locations_idx.push_back(i);
  }


  constructed_trajectory res = construction_heuristic(
    scheduled_locations_idx,
    unscheduled_locations_idx,
    env_state_config);

  std::vector<int> sched_loc = std::get<3>(res);
  MultiWaypointTrajectory mvt = std::get<0>(res);



  for (int i = 1; i < sched_loc.size()-1; i++) {
//    std::cout  << std::endl;
//    std::vector<int> to_try{sched_loc[i-1], sched_loc[i], sched_loc[i+1]};
//    for (auto lalas : to_try) {
//      std::cout << "to_try -> " << lalas << std::endl;
//    }
//    std::cout < << std::endl;
//    auto to_try_optimal = calculate_trajectory_cost_and_optimal_velocities(to_try,
//                                                                 env_state_config);
//    std::cout << "Time now -> " << mvt[i-1].time() + mvt[i].time() << " | time Optimal -> " << std::get<1>(to_try_optimal) << std::endl;
  }


}

int main(int argc, char** argv) {
//  test_pmm(argc, argv);
//  std::string config_file = "config.yaml";
  get_positions_travel_costs("new_config.yaml");
  return 0;
}
