#pragma once

#include "env_config.hpp"
#include "misc_helpers.hpp"
#include <unordered_map>
#include <algorithm>
#include <random>
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

constructed_trajectory run_improved_heuristic(const EnvConfig& env_state_config,
                                           int random_seed,
                                           Scalar cost_leeway_coeff = 1);

constructed_trajectory run_basic_paper_heuristic(const EnvConfig& env_state_config,
                                           int random_seed,
                                           bool optimize_at_every_improvement = false);

constructed_trajectory construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  const EnvConfig& env_params,
  Scalar cost_leeway_coeff = 1,
  MultiWaypointTrajectory mwp_trajectory = MultiWaypointTrajectory{});

constructed_trajectory construction_heuristic_full_optimization(
  std::vector<int> scheduled_locations_idx,
  const EnvConfig& env_params);

// Combination of DH1 and DH2. Remove location for which
// "ratio_of_insertion / max(abs(optimal_flight_time - current_flight_time))" is the lowest.
void destruction_heuristic_3(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             std::vector<Scalar> ratios,
                             const EnvConfig& env_state_config);

// Remove location with most non-optimal connection of heading angle and velocity.
void destruction_heuristic_2(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             const EnvConfig& env_state_config);

// Remove location with lowest insertion ratio.
void destruction_heuristic_1(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             std::vector<Scalar>& ratios);


std::vector<Scalar> calculate_heuristic_ratio(std::vector<int>& scheduled_locations_idx,
                                              MultiWaypointTrajectory& current_trajectory,
                                              std::vector<Scalar> rewards,
                                              const travel_cost_map& precalculated_costs);

std::vector<int> destruction_heuristic_paper(constructed_trajectory& constr_tr,
                                             Scalar percentage,
                                             const EnvConfig& env_params,
                                             std::mt19937_64& rng);


} // namespace agi