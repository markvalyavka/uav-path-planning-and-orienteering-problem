#include "base_heuristic.hpp"

namespace agi {

int randint(int Min, int Max) {
  return std::rand() % (Max + 1 - Min) + Min;
}

constructed_trajectory run_paper_heuristic(EnvConfig& env_state_config,
                                           Scalar cost_leeway_coeff) {

  // idx == location by idx in `location_positions`
  std::vector<int> scheduled_locations_idx = {0, (int)env_state_config.location_positions.size()-1};

  constructed_trajectory initial_constr = construction_heuristic(
    scheduled_locations_idx,
    env_state_config);
  constructed_trajectory best_constr_yet = initial_constr;
  MultiWaypointTrajectory best_tr_yet = std::get<0>(initial_constr);
  Scalar best_cost = std::get<1>(initial_constr);
  Scalar best_reward_yet = std::get<2>(initial_constr);
  std::vector<int> best_scheduled_positions = std::get<3>(initial_constr);


  for (int j = 0; j < 100; j++) {
    std::cout << "First Construction (50%) #" << j << std::endl;

    scheduled_locations_idx = destruction_heuristic_paper(initial_constr, 50, env_state_config);
    initial_constr = construction_heuristic(scheduled_locations_idx, env_state_config, cost_leeway_coeff);
//    break;
    if (std::get<2>(initial_constr) > best_reward_yet) {
      best_constr_yet = initial_constr;
      best_tr_yet = std::get<0>(initial_constr);
      best_cost = std::get<1>(initial_constr);
      best_reward_yet = std::get<2>(initial_constr);
      best_scheduled_positions = std::get<3>(initial_constr);
    }
  }

//  for (int j = 0; j < 100; j++) {
//    std::cout << "Second Construction (20%) #" << j << std::endl;
//
//    scheduled_locations_idx = destruction_heuristic_paper(initial_constr, 20, env_state_config);;
//    initial_constr = construction_heuristic(scheduled_locations_idx, env_state_config);
//    if (std::get<2>(initial_constr) > best_reward_yet) {
//      best_constr_yet = initial_constr;
//      best_tr_yet = std::get<0>(initial_constr);
//      best_cost = std::get<1>(initial_constr);
//      best_reward_yet = std::get<2>(initial_constr);
//      best_scheduled_positions = std::get<3>(initial_constr);
//    }
//  }

  std::cout << "------------------ HERE TEST RESULT  ---------------" << std::endl;
  std::cout << "Best cost -> " << best_cost << std::endl;
  std::cout << "Best reward -> " << best_reward_yet << std::endl;
//  exit(1);
  return best_constr_yet;
}

constructed_trajectory construction_heuristic(
  std::vector<int> scheduled_locations_idx,
  EnvConfig& env_params,
  Scalar cost_leeway_coeff,
  MultiWaypointTrajectory mwp_trajectory) {

  // location_positions  &
  // rewards  &
  // t_max  &
  std::vector<Vector<3>>& location_positions = env_params.location_positions;
  std::vector<Scalar>& rewards = env_params.rewards;
  Scalar t_max = env_params.t_max;

  // precalculated_costs  &
  // velocity_norm_samples  &
  // heading_angle_samples  &
  travel_cost_map& precalculated_costs = env_params.precalculated_costs;
  std::vector<Scalar>& velocity_norm_samples = env_params.velocity_norm_samples;
  std::vector<Scalar>& heading_angle_samples = env_params.heading_angle_samples;


  std::vector<Vector<3>> scheduled_locations{};
  for (int idx : scheduled_locations_idx) {
    scheduled_locations.push_back(env_params.location_positions[idx]);
  }
  std::vector<int> unscheduled_locations_idx = get_missing_values_in_range(scheduled_locations_idx,
                                                                           scheduled_locations_idx[0],
                                                                           scheduled_locations_idx[scheduled_locations_idx.size()-1]);

  MultiWaypointTrajectory current_trajectory;
  if (mwp_trajectory.size() > 0) {
    current_trajectory = mwp_trajectory;
    std::cout << "[NOTE]: trajectory received as param. " << std::endl;
  } else {
    auto current_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
      scheduled_locations_idx,
      env_params,
      true);
    current_trajectory = std::get<0>(current_traj_and_time);
  }


  Scalar current_cost = get_mwp_trajectory_cost(current_trajectory);
  Scalar collected_reward = get_mwp_trajectory_reward(scheduled_locations_idx, rewards);

  while (current_cost < t_max) {

    // <what_idx_to_insert, where_to_insert, velocity>
    std::tuple<int, int, MultiWaypointTrajectory, Scalar, std::vector<Vector<3>>> best_insertion_so_far{};
    Scalar ratio_of_best_insertion_so_far = -1;
    // Try to schedule every unscheduled location.
    for (int unscheduled_idx : unscheduled_locations_idx) {
      //      std::cout << unscheduled_idx << std::endl;
      // For every possible insertion_idx
      //                     possible insertion spot
      //                          v
      // Example: {start_location, end_location}
      for (int insertion_idx = 1; insertion_idx < scheduled_locations.size(); insertion_idx++) {
        //        std::cout << "insertion_idx -> " << insertion_idx << std::endl;
        Scalar curr_min_insertion_cost = DBL_MAX;
        int pred_idx = scheduled_locations_idx[insertion_idx-1];
        int succ_idx = scheduled_locations_idx[insertion_idx];

        // For every possible combination of norm and heading_angle
        for (Scalar norm1 : velocity_norm_samples) {
          for (Scalar angle1 : heading_angle_samples) {
            Vector<3> position_to_schedule = location_positions[unscheduled_idx];
            Scalar pred_to_curr_cost, curr_to_succ_cost, pred_to_succ_cost;

            // Predecessor -> Current [cost]
            Scalar pred_norm = current_trajectory[insertion_idx-1].inp_from_v_norm;
            Scalar pred_angle = current_trajectory[insertion_idx-1].inp_from_v_angle;
            pred_to_curr_cost = precalculated_costs[pred_idx][unscheduled_idx][pred_norm][norm1][pred_angle][angle1];
            if (pred_to_curr_cost == 0) {
              QuadState pred_state;
              QuadState curr_state;
              pred_state.setZero();
              curr_state.setZero();
              pred_state.p = location_positions[pred_idx];
              curr_state.p = position_to_schedule;
              pred_state.v = current_trajectory[insertion_idx-1].get_start_state().v;
              curr_state.v = to_velocity_vector(norm1, angle1);
              PointMassTrajectory3D tr(pred_state, curr_state, env_params.max_acc_per_axis, true);
              if (!tr.exists()) {
//                std::cout << "Not-existing here -> time -> " << tr.time() << " speed -> " << pred_state.v.transpose() << std::endl;
                pred_to_curr_cost = MAX_SCALAR;
              } else {
                pred_to_curr_cost = tr.time();
              }
            }

            // Current -> Successor [cost]
            Scalar succ_norm = current_trajectory[insertion_idx-1].inp_to_v_norm;
            Scalar succ_angle = current_trajectory[insertion_idx-1].inp_to_v_angle;
            curr_to_succ_cost = precalculated_costs[unscheduled_idx][succ_idx][norm1][succ_norm][angle1][succ_angle];
            if (curr_to_succ_cost == 0) {
              QuadState curr_state;
              QuadState succ_state;
              curr_state.setZero();
              succ_state.setZero();
              curr_state.p = position_to_schedule;
              succ_state.p = location_positions[succ_idx];
              curr_state.v = to_velocity_vector(norm1, angle1);
              succ_state.v = current_trajectory[insertion_idx-1].get_end_state().v;
              PointMassTrajectory3D tr(curr_state, succ_state, env_params.max_acc_per_axis, true);
              if (!tr.exists()) {
//                std::cout << "Not-existing here -> time -> " << tr.time() << " speed -> " << succ_state.v.transpose() << std::endl;
                curr_to_succ_cost = MAX_SCALAR;
              } else {
                curr_to_succ_cost = tr.time();
              }
            }

            // Predecessor -> Successor [cost]
            pred_to_succ_cost = current_trajectory[insertion_idx-1].time();
            if (pred_to_curr_cost == 0 || curr_to_succ_cost == 0) {
              std::cout << " This shouldn't happen !!!!!!!!" << std::endl;
            }

            Scalar cost_of_insertion = pred_to_curr_cost + curr_to_succ_cost - pred_to_succ_cost;
            Scalar ratio = rewards[unscheduled_idx] / cost_of_insertion;

            if (ratio > ratio_of_best_insertion_so_far) {

              std::vector<Vector<3>> potential_scheduled_locations = scheduled_locations;
              std::vector<int> potential_scheduled_locations_idx = scheduled_locations_idx;

              potential_scheduled_locations.insert(potential_scheduled_locations.begin()+insertion_idx, position_to_schedule);
              potential_scheduled_locations_idx.insert(potential_scheduled_locations_idx.begin()+insertion_idx, unscheduled_idx);
//              auto try_traj_and_time = calculate_trajectory_cost_and_optimal_velocities(
//                try_scheduled_locations_idx,
//                env_params,
//                true);
//              auto new_cost = std::get<1>(try_traj_and_time);
              Scalar new_cost = cost_of_insertion + current_cost;

              if (new_cost < t_max * cost_leeway_coeff) {

                auto new_trajectory_and_time = calculate_trajectory_cost_and_optimal_velocities(
                  potential_scheduled_locations_idx,
                  env_params,
                  true);

                MultiWaypointTrajectory new_trajectory = std::get<0>(new_trajectory_and_time);
                ratio_of_best_insertion_so_far = ratio;
//                Vector<3> vel_vector = to_velocity_vector(norm1, angle1);
                best_insertion_so_far = {unscheduled_idx, insertion_idx, new_trajectory, std::get<1>(new_trajectory_and_time), potential_scheduled_locations};
              }
            }
          }
        }
      }
    }
    if (std::get<0>(best_insertion_so_far) == 0 || unscheduled_locations_idx.size() == 0)  {
      // If no new points can be scheduled OR all possible points are scheduled, break.
      break;
    }

    //    std::cout << "Best ration -> " << ratio_of_best_insertion_so_far << std::endl;
    int best_loc_idx_to_schedule = std::get<0>(best_insertion_so_far);
    int best_insertion_spot = std::get<1>(best_insertion_so_far);

    current_trajectory = std::get<2>(best_insertion_so_far);
    current_cost = std::get<3>(best_insertion_so_far);
    scheduled_locations = std::get<4>(best_insertion_so_far);
    collected_reward += rewards[best_loc_idx_to_schedule];

    // Update `scheduled_locations_idx` and Remove scheduled_location from `unscheduled_locations_idx`
    scheduled_locations_idx.insert(scheduled_locations_idx.begin()+best_insertion_spot, best_loc_idx_to_schedule);
    unscheduled_locations_idx.erase(std::remove(unscheduled_locations_idx.begin(), unscheduled_locations_idx.end(), best_loc_idx_to_schedule), unscheduled_locations_idx.end());
  }

  std::cout <<  "-------------------------" << std::endl;
  std::cout << "Final cost: " << current_cost << std::endl;
  std::cout << "Collected reward: " << collected_reward << std::endl;
  std::cout << "Actual cost: " << get_mwp_trajectory_cost(current_trajectory) << std::endl;
  std::cout << "Actuual reward: " << get_mwp_trajectory_reward(scheduled_locations_idx, rewards) << std::endl;

  // OUTPUT: MultiWaypointTrajectory, trajectory_cost, trajectory_reward, scheduled_locations_idx, unscheduled_locations_idx
  return {current_trajectory, current_cost, collected_reward, scheduled_locations_idx, unscheduled_locations_idx};
}

void destruction_heuristic_3(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             std::vector<Scalar> ratios,
                             EnvConfig& env_state_config) {

  int min_idx = 1;
  Scalar min_heu = MAX_SCALAR;
  for (int i = 1; i < sched_loc.size()-1; i++) {

    std::vector<int> to_try{sched_loc[i-1], sched_loc[i], sched_loc[i+1]};
    auto to_try_optimal = calculate_trajectory_cost_and_optimal_velocities(to_try,
                                                                           env_state_config,
                                                                           true);
    Scalar optimal_time = std::get<1>(to_try_optimal);
    Scalar curr_time = curr_traj[i-1].time() + curr_traj[i].time();
    Scalar diff = abs(optimal_time - curr_time);
    Scalar heu = ratios[i] / diff;
    if (heu < min_heu) {
      min_heu = heu;
      min_idx = i;
    }
  }
  unsched_loc.push_back(sched_loc[min_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[min_idx]), sched_loc.end());
}

void destruction_heuristic_2(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             MultiWaypointTrajectory& curr_traj,
                             EnvConfig& env_state_config) {

  int max_idx = 1;
  Scalar max_diff = -1;

  for (int i = 1; i < sched_loc.size()-1; i++) {
    std::vector<int> to_try{sched_loc[i-1], sched_loc[i], sched_loc[i+1]};
    auto to_try_optimal = calculate_trajectory_cost_and_optimal_velocities(to_try,
                                                                           env_state_config,
                                                                           true);
    Scalar optimal_time = std::get<1>(to_try_optimal);
    Scalar curr_time = curr_traj[i-1].time() + curr_traj[i].time();
    Scalar diff = abs(optimal_time - curr_time);
    if (diff < max_diff) {
      max_diff = diff;
      max_idx = i;
    }
  }
  unsched_loc.push_back(sched_loc[max_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[max_idx]), sched_loc.end());
}


void destruction_heuristic_1(std::vector<int>& sched_loc,
                             std::vector<int>& unsched_loc,
                             std::vector<Scalar>& ratios) {

  int min_idx = 1;
  Scalar curr_min_ratio = ratios[min_idx];
  for (int i = 1; i < sched_loc.size()-1; i++) {
    if (ratios[i] < curr_min_ratio) {
      curr_min_ratio = ratios[i];
      min_idx = i;
    }
  }
  //  std::cout << "To remove loc -> " << sched_loc[min_idx] << " with ration of " << ratios[min_idx] << std::endl;
  unsched_loc.push_back(sched_loc[min_idx]);
  sched_loc.erase(std::remove(sched_loc.begin(), sched_loc.end(), sched_loc[min_idx]), sched_loc.end());

}

std::vector<Scalar> calculate_heuristic_ratio(std::vector<int>& scheduled_locations_idx,
                                              MultiWaypointTrajectory& current_trajectory,
                                              std::vector<Scalar> rewards,
                                              travel_cost_map& precalculated_costs) {

  std::vector<Scalar> final_ratios{-1};
  // Finding final ratios
  for (int i = 1; i < scheduled_locations_idx.size()-1; i++) {
    Scalar pred_to_curr_cost, curr_to_succ_cost, pred_to_succ_cost;
    Scalar pred_norm = current_trajectory[i-1].inp_from_v_norm;
    Scalar pred_angle = current_trajectory[i-1].inp_from_v_angle;
    Scalar curr_norm = current_trajectory[i-1].inp_to_v_norm;
    Scalar curr_angle = current_trajectory[i-1].inp_to_v_angle;
    Scalar succ_norm = current_trajectory[i].inp_from_v_norm;
    Scalar succ_angle = current_trajectory[i].inp_from_v_angle;

    // Predecessor -> Current [cost]
    pred_to_curr_cost = precalculated_costs[scheduled_locations_idx[i-1]][scheduled_locations_idx[i]][pred_norm][curr_norm][pred_angle][curr_angle];
    // Current -> Successor [cost]
    curr_to_succ_cost = precalculated_costs[scheduled_locations_idx[i]][scheduled_locations_idx[i+1]][curr_norm][succ_norm][curr_angle][succ_angle];
    // Predecessor -> Successor [cost]
    pred_to_succ_cost = precalculated_costs[scheduled_locations_idx[i-1]][scheduled_locations_idx[i+1]][pred_norm][succ_norm][pred_angle][succ_angle];

    Scalar cost_of_insertion = pred_to_curr_cost + curr_to_succ_cost - pred_to_succ_cost;
    Scalar ratio = rewards[scheduled_locations_idx[i]] / cost_of_insertion;
    final_ratios.push_back(ratio);
  }
  final_ratios.push_back(-1);

  return final_ratios;
}


std::vector<int> destruction_heuristic_paper(constructed_trajectory& constr_tr,
                                             Scalar percentage,
                                             EnvConfig& env_params) {
  // ENV PARAMS
  travel_cost_map& travel_costs = env_params.precalculated_costs;
  std::vector<Scalar>& rewards = env_params.rewards;

  MultiWaypointTrajectory mvt = std::get<0>(constr_tr);
  Scalar cost = std::get<1>(constr_tr);
  Scalar reward = std::get<2>(constr_tr);
  std::vector<int> sched_loc = std::get<3>(constr_tr);
  std::vector<int> unsched_loc = std::get<4>(constr_tr);
  std::vector<Scalar> ratios = calculate_heuristic_ratio(sched_loc, mvt, rewards, travel_costs);

  int num_positions_to_remove = sched_loc.size() * percentage / 100;

  //  std::cout << "Percentage -> " << percentage <<"% " << "removes " << num_positions_to_remove << " positions" << std::endl;

  for (int i = 0; i < num_positions_to_remove; i++) {
    int random_number = randint(1, 3);
    if (random_number == 1) {
      destruction_heuristic_1(sched_loc, unsched_loc, ratios);
    } else if (random_number == 2) {
      destruction_heuristic_2(sched_loc, unsched_loc, mvt, env_params);
    } else if (random_number == 3) {
      destruction_heuristic_3(sched_loc, unsched_loc, mvt, ratios, env_params);
    }
    auto new_traj_time = calculate_trajectory_cost_and_optimal_velocities(sched_loc,env_params, true);
    mvt = std::get<0>(new_traj_time);
    cost = std::get<1>(new_traj_time);
    //    std::cout << cost << std::endl;
    ratios = calculate_heuristic_ratio(sched_loc, mvt, rewards, travel_costs);
  }

  return sched_loc;
}

} // namespace agi