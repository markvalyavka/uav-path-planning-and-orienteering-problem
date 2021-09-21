
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <csignal>
#include <fstream>
#include <iostream>
#include <string>

#include "pmm_trajectory3d.hpp"
#include "velocity_search_graph.hpp"
#include "yaml-cpp/yaml.h"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received" << std::endl;
  exit(sig);
}


int test_pmm(int argc, char **argv) {
  // register singal for killing
  std::signal(SIGINT, signal_callback);
  std::cout << "Testing PMM trajectories " << std::endl;

  std::string config_file = "config.yaml";
  YAML::Node config = YAML::LoadFile(config_file);

  QuadState from;
  from.setZero();
  from.p = Vector<3>(0, 0, 1);
  QuadState to;
  to.setZero();
  to.p = Vector<3>(5, 5, 1);
  PointMassTrajectory3D test(from, to, 25);
  std::cout << test << std::endl;

  // -+ value of initial samples
  const Scalar max_yaw_pitch_ang = 15.0;
  // distance between initial samples
  const Scalar precision_yaw_pitch_ang = 15.0;
  // maximal values of the samples
  const Scalar yaw_pitch_cone_angle_boundary = 60.0;

  // velocity norm samples
  const Scalar min_velocity_norm = 1.0;
  const Scalar min_velocity_norm_boundary = 1.0;
  const Scalar max_velocity_norm = 17.0;
  const Scalar precision_velocity_norm = 8.0;

  const Scalar max_acc = 25;

  VelocitySearchGraph vel_search_graph(
    max_acc, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  Vector<3> start_velocity;
  Vector<3> end_velocity;
  Vector<3> start_position;
  Vector<3> end_position;
  config["start"]["velocity"] >> start_velocity;
  config["end"]["velocity"] >> start_velocity;
  config["start"]["position"] >> start_position;
  config["end"]["position"] >> end_position;
  const bool end_free = true;
  std::vector<Vector<3>> gates_waypoints;
  std::vector<Scalar> gates_yaw_deg;

  if (!parseArrayParam<Vector<3>>(config, "gates", gates_waypoints))
    std::cerr << "can not load param gates" << std::endl;

  if (!parseArrayParam<Scalar>(config, "gates_orientations", gates_yaw_deg))
    std::cerr << "can not load param gates_orientations" << std::endl;

  gates_waypoints.insert(gates_waypoints.begin(), start_position);
  gates_waypoints.push_back(end_position);

  std::cout << "gates_waypoints.size() " << gates_waypoints.size() << std::endl;
  std::cout << "gates_yaw_deg.size() " << gates_yaw_deg.size() << std::endl;


  gates_waypoints.resize(3);
  gates_yaw_deg.resize(3);
  vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, end_velocity, gates_yaw_deg, end_free);

  // save_track_trajectory(trajectories, shortest_time,
  //                       output_folder_ + "samples_pmm.csv");

  // save_track_trajectory_equidistant(single, trajectories[0].time(),
  //                                   "samples_equidistant.csv");

  return 0;
}

int main(int argc, char **argv) {
  test_pmm(argc, argv);
  return 0;
}
