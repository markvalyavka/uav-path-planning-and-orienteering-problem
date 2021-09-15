
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <csignal>
#include <fstream>
#include <iostream>
#include <string>

#include "pmm_trajectory3d.hpp"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received" << std::endl;
  exit(sig);
}

int test_pmm(int argc, char **argv) {
  // register singal for killing
  std::signal(SIGINT, signal_callback);
  std::cout << "Testing PMM trajectories " << std::endl;
  QuadState from;
  from.setZero();
  from.p = Vector<3>(0, 0, 1);
  QuadState to;
  to.setZero();
  to.p = Vector<3>(5, 5, 1);
  PointMassTrajectory3D test(from, to, 25);
  std::cout << test << std::endl;


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
