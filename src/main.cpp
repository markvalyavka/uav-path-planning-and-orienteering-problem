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
#include "base_heuristic.hpp"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received. Exiting.. " << std::endl;
  exit(sig);
}


MultiWaypointTrajectory optimize_with_cone_refocusing(std::vector<Vector<3>> gates_waypoints,
             std::vector<Scalar> gates_vel_norms,
             std::vector<Scalar> gates_yaw_deg,
             Vector<3> start_velocity) {

  // Default cone refocusing params.
  const Scalar min_velocity_norm = 0;
  const Scalar max_velocity_norm = 2.12132;
  const Scalar min_velocity_norm_boundary = 0.01;
  const Scalar precision_velocity_norm = 0.1;  // Distance between velocity norms
  const Scalar max_acc_norm = 1.06066;

  const Scalar yaw_pitch_cone_angle_boundary = 40.0;
  const Scalar max_yaw_pitch_ang = 20.0;
  const Scalar precision_yaw_pitch_ang = 20.0; // Distance between yaw angles


  const bool sample_start_velocity = true;

  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  std::vector<Scalar> gates_pitch_deg(gates_yaw_deg.size(), 0);

  MultiWaypointTrajectory refocused_trajectory = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, gates_yaw_deg,
    gates_pitch_deg, gates_vel_norms, sample_start_velocity);

  return refocused_trajectory;
}

void get_positions_travel_costs(std::string config_file, int argc, char** cli_args)
{
  // LOAD CONFIG
  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();
//  exit(1);

  std::cout << "------------------  1. FINAL STATS AFTER RUN PAPER HEU  ---------------" << std::endl;
  constructed_trajectory heuristic_result = run_paper_heuristic(env_state_config);
  MultiWaypointTrajectory best_tr_yet = std::get<0>(heuristic_result);
  Scalar best_cost = std::get<1>(heuristic_result);
  Scalar best_reward_yet = std::get<2>(heuristic_result);
  std::vector<int> best_scheduled_positions = std::get<3>(heuristic_result);

  // -----------------------------------------------------------------------
  MultiWaypointTrajectory final_trajectory = best_tr_yet;
  Scalar final_cost = best_cost;
  Scalar final_reward = best_reward_yet;
  std::vector<int> final_scheduled_positions_idx = best_scheduled_positions;
  // -----------------------------------------------------------------------

  MultiWaypointTrajectory current_trajectory = best_tr_yet;
  Scalar current_cost = best_cost;
  Scalar current_reward = best_reward_yet;
  std::vector<int> current_scheduled_positions_idx = best_scheduled_positions;
  std::vector<Vector<3>> current_scheduled_positions{};


  int imp_iterations = 6;


  for(int imp_i = 0; imp_i < imp_iterations; imp_i++) {
    std::cout << "------------------ CONE REFOCUSING (Iter #" << imp_i << ") ---------------" << std::endl;
    std::vector<Scalar> vel_norms = get_mwp_trajectory_velocities(current_trajectory);
    std::vector<Scalar> yaw_angles = get_mwp_trajectory_yaw_angles(current_trajectory);
    std::vector<Vector<3>> scheduled_positions;
    for (int idx : current_scheduled_positions_idx) {
      scheduled_positions.push_back(env_state_config.location_positions[idx]);
    }
    Vector<3> start_vel = to_velocity_vector(current_trajectory[0].inp_from_v_norm, current_trajectory[0].inp_from_v_angle);
    MultiWaypointTrajectory cone_refocused_trajectory = optimize_with_cone_refocusing(scheduled_positions, vel_norms, yaw_angles, start_vel);


    current_trajectory = cone_refocused_trajectory;
    current_cost = get_mwp_trajectory_cost(cone_refocused_trajectory);
    current_reward = get_mwp_trajectory_reward(current_scheduled_positions_idx, env_state_config.rewards);
    current_scheduled_positions = scheduled_positions;

    if (current_cost > env_state_config.t_max) {
      // If trajectory cost after refocusing is greater than budget,
      // we can stop the algorithm as there is no chance of further improvement.
      break;
    }
    if (current_reward > final_reward) {
      final_trajectory = current_trajectory;
      final_cost = current_cost;
      final_reward = current_reward;
      final_scheduled_positions_idx = current_scheduled_positions_idx;
    }

    std::cout << "------------------ AUGMENTING (Iter #" << imp_i << ") ---------------" << std::endl;
    constructed_trajectory augmented_trajectory = construction_heuristic(
      current_scheduled_positions_idx,
      env_state_config,
      current_trajectory);
    current_trajectory = std::get<0>(augmented_trajectory);
    current_cost = std::get<1>(augmented_trajectory);
    current_reward = std::get<2>(augmented_trajectory);
    current_scheduled_positions_idx = std::get<3>(augmented_trajectory);
    std::vector<Vector<3>> temp_scheduled_positions{};
    for (int idx : current_scheduled_positions_idx) {
      temp_scheduled_positions.push_back(env_state_config.location_positions[idx]);
    }
    current_scheduled_positions = temp_scheduled_positions;
    if (current_cost < env_state_config.t_max & current_reward > final_reward) {
      final_trajectory = current_trajectory;
      final_cost = current_cost;
      final_reward = current_reward;
      final_scheduled_positions_idx = current_scheduled_positions_idx;
    }
  }

  std::cout << "------------------ FINAL RESULT  ---------------" << std::endl;
  std::cout << "Final cost -> " << final_cost << std::endl;
  std::cout << "Final reward -> " << final_reward << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(final_trajectory, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;
//  exit(1);
//  print_detailed_mwp_stats(final_trajectory, env_state_config.max_acc_per_axis);

}

int main(int argc, char** argv) {

  // 1. Add `leeway` to construction heuristic.
  // 2. Heuristic by variably changing `leeway` (progressively go from 1.10 -> 1.0)?
  // 3. Try to optimize final solution with very powerful (lots of samples cones refocusing?)

  srand(3);
  get_positions_travel_costs("/Users/markv/pmm_planner/new_config.yaml", argc, argv);

  return 0;
}
