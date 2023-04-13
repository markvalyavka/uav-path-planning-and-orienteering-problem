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

#include <thread>
#include <future>

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
  std::vector<Scalar> gates_pitch_deg(gates_yaw_deg.size(), 0);

  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);
  MultiWaypointTrajectory refocused_trajectory = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, gates_yaw_deg,
    gates_pitch_deg, gates_vel_norms, sample_start_velocity);
  return refocused_trajectory;
}


std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>> run_improved_trajectory_algorithm_with_cost_coeff(EnvConfig& env_state_config,
                                                                                                                        int random_seed,
                                                                                                                        Scalar cost_leeway_coeff) {

  std::cout << "------------------  1. FINAL STATS AFTER RUN PAPER HEU  ---------------" << std::endl;
  constructed_trajectory heuristic_result = run_paper_heuristic(env_state_config, random_seed, cost_leeway_coeff);
  // -----------------------------------------------------------------------
  MultiWaypointTrajectory final_trajectory = std::get<0>(heuristic_result);
  Scalar final_cost = std::get<1>(heuristic_result);
  Scalar final_reward = std::get<2>(heuristic_result);
  std::vector<int> final_scheduled_positions_idx = std::get<3>(heuristic_result);
  // -----------------------------------------------------------------------
//  std::cout << "initial cost -> " << final_cost << std::endl;
//  std::cout << "initial reward -> " << final_reward << std::endl;
  MultiWaypointTrajectory current_trajectory = std::get<0>(heuristic_result);
  Scalar current_cost = std::get<1>(heuristic_result);
  Scalar current_reward = std::get<2>(heuristic_result);
  std::vector<int> current_scheduled_positions_idx = std::get<3>(heuristic_result);

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
    if (current_trajectory.size() < 2) {
      break;
    }
    current_cost = get_mwp_trajectory_cost(cone_refocused_trajectory);
    current_reward = get_mwp_trajectory_reward(current_scheduled_positions_idx, env_state_config.rewards);
    if (current_cost > env_state_config.t_max || current_cost <= 0) {
      // If trajectory cost after refocusing is greater than budget,
      // we can stop the algorithm as there is no chance of further improvement.
//      std::cout << "hre" << std::endl;
      break;
    }
    if (current_reward >= final_reward) {
      final_trajectory = current_trajectory;
      final_cost = current_cost;
      final_reward = current_reward;
      final_scheduled_positions_idx = current_scheduled_positions_idx;
    }

    std::cout << "------------------ AUGMENTING (Iter #" << imp_i << ") ---------------" << std::endl;
    constructed_trajectory augmented_trajectory = construction_heuristic(
      current_scheduled_positions_idx,
      env_state_config,
      cost_leeway_coeff,
      current_trajectory);

    current_trajectory = std::get<0>(augmented_trajectory);
    current_cost = std::get<1>(augmented_trajectory);
    current_reward = std::get<2>(augmented_trajectory);
    current_scheduled_positions_idx = std::get<3>(augmented_trajectory);
    if (current_cost < env_state_config.t_max & current_reward > final_reward) {
      final_trajectory = current_trajectory;
      final_cost = current_cost;
      final_reward = current_reward;
      final_scheduled_positions_idx = current_scheduled_positions_idx;
    }
  }

//  std::cout << "IN FUNCT REWARD -> " << final_reward << std::endl;

  return {final_trajectory, final_cost, final_reward, final_scheduled_positions_idx};
}


void run_improved_trajectory_algorithm(std::string config_file, int random_seed) {
  // LOAD CONFIG
  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();

  // --------------------------------------------------------
  MultiWaypointTrajectory final_trajectory{};
  Scalar final_cost = 0;
  Scalar final_reward = 0;
  std::vector<int> final_scheduled_locations_idx{};
  // --------------------------------------------------------

  // Create a vector to store the results.
  std::vector<std::future<std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>>>> future_results;
  std::vector<std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>>> output;


  std::vector<Scalar> cost_leeway_coeffs = {1, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.16, 1.3, 1.5};
//  std::vector<Scalar> cost_leeway_coeffs = {1.02};
  // Create a thread for each value of cost_leeway_coeff.
  for (const auto cost_leeway_coeff : cost_leeway_coeffs) {
    future_results.emplace_back(std::async(std::launch::async,
                                           run_improved_trajectory_algorithm_with_cost_coeff,
                                           std::ref(env_state_config), random_seed, cost_leeway_coeff));
  }
  for (auto& fut_res : future_results) {
    output.push_back(fut_res.get());
  }

  int i = 0;
  for (auto& result : output) {
    Scalar curr_cost = std::get<1>(result);
    Scalar curr_reward = std::get<2>(result);
    std::cout << "coeff -> " << cost_leeway_coeffs[i] << " cost -> " << curr_cost << " reward -> " << curr_reward << std::endl;
    if (curr_cost < env_state_config.t_max & curr_reward > final_reward) {
      final_trajectory = std::get<0>(result);
      final_cost = std::get<1>(result);
      final_reward = std::get<2>(result);
      final_scheduled_locations_idx = std::get<3>(result);
    }
    i++;
  }

  std::cout << "------------------ FINAL RESULT [IMPROVED_ALGORITHM]  ---------------" << std::endl;
  std::cout << "Final cost -> " << final_cost << std::endl;
  std::cout << "Final reward -> " << final_reward << std::endl;
  save_trajectory_results(final_trajectory, final_reward);
  //  exit(1);
  //  print_detailed_mwp_stats(final_trajectory, env_state_config.max_acc_per_axis);
}

void run_basic_trajectory_algorithm(std::string config_file, int random_seed) {
  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();


  constructed_trajectory heuristic_result = run_paper_heuristic(env_state_config, random_seed, 1);
  // -----------------------------------------------------------------------
  MultiWaypointTrajectory final_trajectory = std::get<0>(heuristic_result);
  Scalar final_cost = std::get<1>(heuristic_result);
  Scalar final_reward = std::get<2>(heuristic_result);
  std::vector<int> final_scheduled_positions_idx = std::get<3>(heuristic_result);
  // -----------------------------------------------------------------------


  std::cout << "------------------ FINAL RESULT [BASIC_ALGORITHM]  ---------------" << std::endl;
  std::cout << "Final cost -> " << final_cost << std::endl;
  std::cout << "Final reward -> " << final_reward << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(final_trajectory, "samples_pmm.csv");
  std::cout << "Saved equitemporal." << std::endl;
}

int main(int argc, char** argv) {


  // 1. Add `leeway` to construction heuristic. (`check`)
  // 2. Heuristic by variably changing `leeway` (progressively go from 1.10 -> 1.0) -> threads? (`check`)
  //   2.1. For small t_max, even 1.1 coeff might be too small (`check`)
  // 3. Try to optimize final solution with very powerful (lots of samples cones refocusing?)
  // 4. Allow `start_vel` sampling in cone_refocus. (from stash)

  int random_seed = 4;

  run_improved_trajectory_algorithm("/Users/markv/pmm_planner/input_configs/cfg_ts3.yaml", random_seed);
//  run_basic_trajectory_algorithm("/Users/markv/pmm_planner/cfg_ts1.yaml", random_seed);

  return 0;
}
