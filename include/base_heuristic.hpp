#pragma once

#include "env_config.hpp"
#include "misc_helpers.hpp"
#include <unordered_map>
#include <algorithm>
#include "gravity.hpp"
#include "pmm_trajectory3d.hpp"
#include "three_acc.hpp"
#include "timer.hpp"
#include "tuples_hash.hpp"
#include "velocity_search_graph.hpp"
#include "yaml-cpp/yaml.h"
#include "trajectory_calculation.hpp"

namespace agi {

typedef std::tuple<MultiWaypointTrajectory, Scalar, Scalar, std::vector<int>, std::vector<int>> constructed_trajectory;

int randint(int Min, int Max);

constructed_trajectory run_paper_heuristic(EnvConfig& env_state_config);

constructed_trajectory construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  EnvConfig& env_params,
  MultiWaypointTrajectory mwp_trajectory = MultiWaypointTrajectory{});

void destruction_heuristic_3(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             std::vector<Scalar> ratios,
                             EnvConfig& env_state_config);


void destruction_heuristic_2(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             EnvConfig& env_state_config);

void destruction_heuristic_1(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             std::vector<Scalar>& ratios);


std::vector<Scalar> calculate_heuristic_ratio(std::vector<int>& scheduled_locations_idx,
                                              MultiWaypointTrajectory& current_trajectory,
                                              std::vector<Scalar> rewards,
                                              travel_cost_map& precalculated_costs);

std::vector<int> destruction_heuristic_paper(const constructed_trajectory& constr_tr,
                                             Scalar percentage,
                                             EnvConfig& env_params);


} // namespace agi