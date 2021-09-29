
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


int test_pmm(int argc, char** argv) {
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

  const Scalar max_acc_norm = 32.94;

  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
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
  Scalar sum_times = 0;
  MultiWaypointTrajectory tr = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, end_velocity, gates_yaw_deg, end_free,
    false);
  std::cout << "output tr size " << tr.size() << std::endl;
  for (size_t i = 0; i < tr.size(); i++) {
    std::cout << i << " vel " << tr[i].get_end_state().v.transpose()
              << std::endl;
    // std::cout << tr[i] << std::endl;
    sum_times += tr[i].time();
    if (tr[i].time() - tr[i].time_min() > PRECISION_PMM_VALUES) {
      std::cout << "bad time!!!!!!" << std::endl;
    }
  }
  std::cout << "total time " << sum_times << std::endl;
  std::vector<QuadState> samples =
    VelocitySearchGraph::getTrajectoryEquidistantStates(tr, 0.8);

  // A*x = b
  // x coefficients
  // A polynomial coefficients
  // b = [p0,p1,v0,v1,a0,a1]
  std::vector<std::vector<Vector<6>>> polynomials;
  std::cout << "samples.size() " << samples.size() << std::endl;
  polynomials.resize(samples.size() - 1);
  for (size_t i = 1; i < samples.size(); i++) {
    QuadState& from = samples[i - 1];
    QuadState& to = samples[i];
    Scalar dt = to.t - from.t;
    // x = [a5,a4,a3,a2,a1,a0]
    for (size_t axi = 0; axi < 3; axi++) {
      Vector<8> b;
      b << from.p(axi), to.p(axi), from.v(axi), to.v(axi), from.a(axi),
        to.a(axi), 0, 0;
      Matrix<8, 6> A;
      Vector<6> tau = Vector<6>::Ones();
      for (int i = 1; i < 6; i++) {
        tau(i) = tau(i - 1) * dt;
      }
      std::cout << "tau " << tau.transpose() << std::endl;
      A.row(0) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;                    // p0
      A.row(1) << tau(5), tau(4), tau(3), tau(2), tau(1), tau(0);  // p1
      A.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;                    // v0
      A.row(3) << 5.0 * tau(4), 4.0 * tau(3), 3.0 * tau(2), 2.0 * tau(1),
        1.0 * tau(0), 0.0;                       // v1
      A.row(4) << 0.0, 0.0, 0.0, 2.0, 0.0, 0.0;  // a0
      A.row(5) << 20.0 * tau(3), 12.0 * tau(2), 6.0 * tau(1), 2.0 * tau(0), 0.0,
        0.0;                                     // a1
      A.row(6) << 0.0, 0.0, 6.0, 0.0, 0.0, 0.0;  // j0
      A.row(7) << 60.0 * tau(2), 24.0 * tau(1), 6.0 * tau(0), 0, 0.0,
        0.0;  // j1
      // A.row(6) << 0.0, 24.0, 0.0, 0.0, 0.0, 0.0;  // j0
      // A.row(7) << 120.0 * tau(1), 24.0 * tau(0), 0.0, 0.0, 0.0,
      //   0.0;  // j1

      std::cout << "solve" << std::endl;
      Vector<6> p = A.colPivHouseholderQr().solve(b);
      polynomials[i - 1].push_back(p);
      std::cout << "solved" << std::endl;
      std::cout << "p" << i << "[" << axi << "] " << p.transpose() << std::endl;
    }
  }

  std::ofstream myfile;
  myfile.open("polynomials.csv");
  if (myfile.is_open()) {
    if (samples.size() > 0 && samples[0].size() > 0) {
      myfile << "i,axi,tfrom,tto,a5,a4,a3,a2,a1,a0" << std::endl;
    }

    for (int var1 = 0; var1 < polynomials.size(); ++var1) {
      const Scalar tfrom = samples[var1].t;
      const Scalar tto = samples[var1 + 1].t;
      for (size_t axi = 0; axi < 3; axi++) {
        const Vector<6>& p = polynomials[var1][axi];
        myfile << var1 << "," << axi << "," << tfrom << "," << tto << ","
               << p(0) << "," << p(1) << "," << p(2) << "," << p(3) << ","
               << p(4) << "," << p(5);
        myfile << std::endl;
      }
    }

    myfile.close();
  }

  std::cout << "out" << std::endl;

  VelocitySearchGraph::saveTrajectoryEquitemporal(tr, "samples_pmm.csv");

  std::cout << "saved equitemporal" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(tr, "samples_equidistant.csv");
  std::cout << "saved equidistant" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(
    tr, "samples_equidistant_08.csv", 0.8);
  std::cout << "saved equidistant 0.8" << std::endl;

  return 0;
}

int main(int argc, char** argv) {
  test_pmm(argc, argv);
  return 0;
}
