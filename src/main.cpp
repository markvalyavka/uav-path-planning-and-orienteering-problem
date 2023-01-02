
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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
#include "yaml-cpp/yaml.h"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received. Exiting.. " << std::endl;
  exit(sig);
}

struct node {
  QuadState state;
  Scalar g = 0;  // current cost
  Scalar h = 0;  // heuristic
  int sample_index;
  node* parent;
  Scalar f() { return g + h; };
};

typedef std::tuple<short, short, short, short, short, short, short, short,
                   short>
  dist_state;

dist_state make_diststate(QuadState state) {
  return std::make_tuple((short)(state.p(0) * 5.0), (short)(state.p(1) * 5.0),
                         (short)(state.p(2) * 5.0), (short)(state.v(0) * 2.0),
                         (short)(state.v(1) * 2.0), (short)(state.v(2) * 2.0),
                         (short)(state.a(0) / 3.0), (short)(state.a(1) / 3.0),
                         (short)(state.a(2) / 3.0));
}


// point mass model
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
  if (!parseArrayParam<Vector<3>>(config, "gates", gates_waypoints))
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

  std::cout << "I got here!" << std::endl;
  // Program stops here

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

std::vector<Scalar> get_velocity_norm_samples(int samples_num, Scalar max_velocity) {
  std::vector<Scalar> samples{};
  for (int i = 1; i < samples_num+1; i++) {
    samples.emplace_back((i-1)*max_velocity/(samples_num-1));
  }
  return samples;
}

std::vector<Scalar> get_heading_angle_samples(int samples_num) {
  std::vector<Scalar> samples{};
  for (int i = 1; i < samples_num+1; i++) {
    samples.emplace_back(2*i*M_PI_2/samples_num);
  }
  return samples;
}

Vector<3> to_velocity_vector(Scalar velocity_norm, Scalar angle) {
  return Vector<3>(velocity_norm*cos(angle), velocity_norm*sin(angle), 0);
}

std::tuple<Scalar, Scalar> to_velocity_norm_and_angle(Vector<3> velocity_vector) {
  Scalar x_squared = velocity_vector[0] * velocity_vector[0];
  Scalar y_squared = velocity_vector[1] * velocity_vector[1];
  Scalar velocity_norm = sqrt(x_squared + y_squared);
  Scalar heading_angle = atan2(velocity_vector[1], velocity_vector[0]);
  return std::make_tuple(velocity_norm, heading_angle);
}


void get_positions_travel_costs()
{
  const Scalar single_axis = 100.0;
  Vector<3> max_acc_per_axis = Vector<3>::Constant(single_axis);
  Scalar max_velocity = 20;
  int V = 8;
  int H = 8;
  std::vector<Scalar> velocity_samples = get_velocity_norm_samples(V, max_velocity);
  std::vector<Scalar> heading_angle_samples = get_heading_angle_samples(H);
  std::cout << "Velocity vector -> " << to_velocity_vector(18, 2.35619) << std::endl;
  Scalar norm, hangle;
  std::tie(norm, hangle) = to_velocity_norm_and_angle(to_velocity_vector(18, 2.35619));
  std::cout << "norm -> " << norm << " angle -> " << hangle << std::endl;
//  for (auto sample: heading_angle_samples) {
//    std::cout << "Heading angle sample -> " << sample << std::endl;
//  }

  QuadState loc_1;
  loc_1.setZero();
//  loc_1.p = Vector<3>(8.34, 6.34, 0.757);
//  loc_1.v = Vector<3>(12.4, 4.53, -2.59);
  loc_1.p = Vector<3>(0, 0, 0);
  loc_1.v = Vector<3>(12, 0, 0);

  QuadState loc_2;
  loc_2.setZero();
//  loc_2.p = Vector<3>(9.09, 6.26, 1.08);
//  loc_2.v = Vector<3>(15.9748, 0, -5.81434);
  loc_2.p = Vector<3>(10, 10, 0);
  loc_2.v = Vector<3>(12, 12, 0);

  QuadState loc_3;
  loc_3.setZero();
//    to.p = Vector<3>(9.09, 6.26, 1.08);
//    to.v = Vector<3>(15.9748, 0, -5.81434);
  loc_3.p = Vector<3>(20, 13, 0);
  loc_3.v = Vector<3>(0, 12, 0);

  QuadState loc_4;
  loc_4.setZero();
  //    to.p = Vector<3>(9.09, 6.26, 1.08);
  //    to.v = Vector<3>(15.9748, 0, -5.81434);
  loc_4.p = Vector<3>(0, 14, 0);
  loc_4.v = Vector<3>(-20, -20, 0);

  PointMassTrajectory3D two_point_trajectory1(loc_1, loc_2, max_acc_per_axis, true);
  PointMassTrajectory3D two_point_trajectory2(loc_2, loc_3, max_acc_per_axis, true);
  PointMassTrajectory3D two_point_trajectory3(loc_3, loc_4, max_acc_per_axis, true);

  //  std::cout << two_point_trajectory.time() << std::endl;
//  std::cout << two_point_trajectory.exists() << std::endl;
  std::vector<PointMassTrajectory3D> traj = {two_point_trajectory1, two_point_trajectory2, two_point_trajectory3};
//  auto time_min = traj[0].time_min();
//  std::cout << "Time min: " << time_min << std::endl;
//  const Scalar traj_length = traj[0].get_length_between_times(0, time_min);

//  std::cout <<  "--------------------" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquitemporal(traj, "samples_pmm.csv");
  std::cout << "saved equitemporal" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(traj,
                                                 "samples_equidistant.csv");
  std::cout << "saved equidistant" << std::endl;
}

int main(int argc, char** argv) {
//  test_pmm(argc, argv);
  std::string config_file = "config.yaml";
  get_positions_travel_costs();
  return 0;
}
