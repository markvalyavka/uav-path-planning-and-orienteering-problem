#pragma once

#include "env_config.hpp"
#include "misc_helpers.hpp"

namespace agi {


std::tuple<MultiWaypointTrajectory, Scalar> calculate_trajectory_cost_and_optimal_velocities(std::vector<int> &scheduled_locations_idx,
                                                                                             EnvConfig& env_params,
                                                                                             bool sample_start_velocity = false);

std::tuple<MultiWaypointTrajectory, Scalar> imp_calculate_trajectory_cost_and_optimal_velocities(std::vector<int> &scheduled_locations_idx,
                                                                                             EnvConfig& env_params,
                                                                                             bool sample_start_velocity = false);


}
