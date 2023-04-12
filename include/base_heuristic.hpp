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

constructed_trajectory run_paper_heuristic(EnvConfig& env_state_config,
                                           Scalar cost_leeway_coeff = 1);

constructed_trajectory construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  EnvConfig& env_params,
  Scalar cost_leeway_coeff = 1,
  MultiWaypointTrajectory mwp_trajectory = MultiWaypointTrajectory{});

// Combination of DH1 and DH2. Remove location for which
// "ratio_of_insertion / max(abs(optimal_flight_time - current_flight_time))" is the lowest.
void destruction_heuristic_3(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             std::vector<Scalar> ratios,
                             EnvConfig& env_state_config);

// Remove location with most non-optimal connection of heading angle and velocity.
void destruction_heuristic_2(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             EnvConfig& env_state_config);

// Remove location with lowest insertion ratio.
void destruction_heuristic_1(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             std::vector<Scalar>& ratios);


std::vector<Scalar> calculate_heuristic_ratio(std::vector<int>& scheduled_locations_idx,
                                              MultiWaypointTrajectory& current_trajectory,
                                              std::vector<Scalar> rewards,
                                              travel_cost_map& precalculated_costs);

std::vector<int> destruction_heuristic_paper(constructed_trajectory& constr_tr,
                                             Scalar percentage,
                                             EnvConfig& env_params);


} // namespace agi