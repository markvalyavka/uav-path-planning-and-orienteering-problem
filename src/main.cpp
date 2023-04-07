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

void get_positions_travel_costs(std::string config_file, int argc, char** cli_args)
{
  // LOAD CONFIG
  EnvConfig env_state_config(config_file);
  env_state_config.generate_samples_with_simple_sampling();
  env_state_config.generate_precalculated_graph_of_costs();


//  std::vector test_scheduled_locs_idx{0, 6, 5, 4, 3, 19, 18, 15, 14, 16, 9, 13, 20};
//
//  auto test_tr_and_time = imp_calculate_trajectory_cost_and_optimal_velocities(test_scheduled_locs_idx,
//                                                   env_state_config,
//                                                                           true);
//  auto test_tr = std::get<0>(test_tr_and_time);
//  auto test_time = std::get<1>(test_tr_and_time);
//  auto test_reward = get_mwp_trajectory_reward(test_scheduled_locs_idx, env_state_config.rewards);
//  Scalar test_actual_time = 0;
//  for (auto tr : test_tr) {
//    test_actual_time += tr.time();
//  }
//  auto test_func_time = get_mwp_trajectory_cost(test_tr);
//  std::cout << "test_actual_time -> " << test_actual_time << std::endl;
//  std::cout << "test_time -> " << test_time << std::endl;
//  std::cout << "test_reward -> " << test_reward << std::endl;
//  std::cout << "test_func_time -> " << test_func_time << std::endl;
//  exit(1);


  std::cout << "------------------  1. FINAL STATS AFTER RUN PAPER HEU  ---------------" << std::endl;
  auto heu_result = run_paper_heuristic(env_state_config);
  MultiWaypointTrajectory best_tr_yet = std::get<0>(heu_result);
  Scalar best_cost = std::get<1>(heu_result);
  Scalar best_reward_yet = std::get<2>(heu_result);
  std::vector<int> best_scheduled_positions = std::get<3>(heu_result);

  std::cout << "Reward: " << best_reward_yet << std::endl;
  std::cout << "Cost: " << best_cost << std::endl;
  std::cout << "Scheduled locations:" << std::endl;
  for (auto sch: best_scheduled_positions) {
    std::cout << sch << " -> ";
  }
  VelocitySearchGraph::saveTrajectoryEquitemporal(best_tr_yet, "samples_pmm.csv");
  std::cout <<  "--------------------" << std::endl;
  std::cout << "Saved equitemporal." << std::endl;


  std::cout << "------------------ 2. START CONE REFOCUSING  ---------------" << std::endl;
  std::vector<Scalar> vel_norms = get_mwp_trajectory_velocities(best_tr_yet);
  std::vector<Scalar> yaw_angles = get_mwp_trajectory_yaw_angles(best_tr_yet);
  std::vector<Vector<3>> scheduled_loc_gates;
  for (int pos_idx : best_scheduled_positions) {
    scheduled_loc_gates.push_back(env_state_config.location_positions[pos_idx]);
  }
  Vector<3> start_vel = to_velocity_vector(best_tr_yet[0].inp_from_v_norm, best_tr_yet[0].inp_from_v_angle);
  auto cone_refoces_tr = test_pmm(scheduled_loc_gates, vel_norms, yaw_angles, start_vel);

  std::cout << "------------------ 3. AUGMENTING  ---------------" << std::endl;
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

  std::cout << "------------------ 4. START CONE REFOCUSING - 2  ---------------" << std::endl;
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
  std::cout << "Saved equitemporal." << std::endl;

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2 22222222222222
  std::cout << "------------------ 5. AUGMENTING  ---------------" << std::endl;

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
//  std::cout << " cost -> " << cost3 << std::endl;
  std::cout << " actual cost -> " << actual_cost3 << std::endl;
//  std::cout << " reward -> " << reward3 << std::endl;
  std::cout << " actual reward -> " << actual_reward3 << std::endl;

  std::cout << "------------------ 6. START CONE REFOCUSING - 3  ---------------" << std::endl;

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
  std::cout << "------------------ 7. FINAL STATS  ---------------" << std::endl;
  std::cout << "Cost -> " << get_mwp_trajectory_cost(cone_refoces_tr3) << std::endl;
  std::cout << "Reward -> " << get_mwp_trajectory_reward(scheduled_idx3, env_state_config.rewards) << std::endl;
  exit(1);
}

int main(int argc, char** argv) {
//  std::cout << to_velocity_vector(1.06066, 1.58) << std::endl;
//  exit(1);
//  test_pmm();
//
  srand(2);
  get_positions_travel_costs("/Users/markv/pmm_planner/new_config.yaml", argc, argv);

  return 0;
}
