
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
  QuadState to;
  PointMassTrajectory3D test(from, to, 25);
  std::cout << test << std::endl;
  return 0;
}

int main(int argc, char **argv) {
  test_pmm(argc, argv);
  return 0;
}
