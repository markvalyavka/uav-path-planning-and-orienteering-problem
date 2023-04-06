#include "misc_helpers.hpp"


namespace agi {


std::vector<Scalar> get_velocity_norm_samples(int samples_num, Scalar max_velocity) {
  std::vector<Scalar> samples{};
  if (samples_num == 1) {
    return std::vector<Scalar>{max_velocity};
  }
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
//  std::cout << "angle -> " << angle << std::endl;
  return Vector<3>(velocity_norm*cos(angle), velocity_norm*sin(angle), 0);
}

std::tuple<Scalar, Scalar> to_velocity_norm_and_angle(Vector<3> velocity_vector) {
  Scalar x_squared = velocity_vector[0] * velocity_vector[0];
  Scalar y_squared = velocity_vector[1] * velocity_vector[1];
  Scalar velocity_norm = sqrt(x_squared + y_squared);
  Scalar heading_angle = atan2(velocity_vector[1], velocity_vector[0]);
  return std::make_tuple(velocity_norm, heading_angle);
}

std::vector<std::tuple<Vector<3>, Scalar, Scalar>> generate_total_sample_tuples(
  norm_angle_to_velocity_vector_map &norm_angle_to_vector,
  std::vector<Scalar> &velocity_norm_samples,
  std::vector<Scalar> &heading_angle_samples) {
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>> samples{};
  for (Scalar norm : velocity_norm_samples) {
    for (Scalar angle : heading_angle_samples) {
      Vector<3> velocity_vector = norm_angle_to_vector[norm][angle];
      samples.push_back({velocity_vector, norm, angle});
    }
  }
  return samples;
}

void populate_norm_angle_to_velocity_vector_map(norm_angle_to_velocity_vector_map &map,
                                                std::vector<Scalar> &vel_norm_samples,
                                                std::vector<Scalar> &angle_samples) {
  for (Scalar norm : vel_norm_samples) {
    for (Scalar angle : angle_samples) {
      Vector<3> velocity_vec = to_velocity_vector(norm, angle);
      //      std::cout << "norm -> " << norm << ", angle -> " << angle
      //            << " = Vector<" << velocity_vec[0] << ", " << velocity_vec[1] << ", " << velocity_vec[2] << "> " << std::endl;
      map[norm][angle] = velocity_vec;
    }
  }
}

void imp_populate_precalculated_travel_costs_map(imp_travel_cost_map& travel_costs,
                                             std::vector<Vector<3>>& location_positions,
                                             std::vector<std::tuple<Vector<3>, Scalar, Scalar>> velocity_samples_tuples,
                                             Vector<3> &max_acc_per_axis) {

  QuadState test_loc_1;
  QuadState test_loc_2;
  test_loc_1.setZero();
  test_loc_2.setZero();
  for (int loc1_id = 0; loc1_id < location_positions.size(); loc1_id++) {
    for (int loc2_id = 0; loc2_id < location_positions.size(); loc2_id++) {
      if (loc1_id == loc2_id) {continue;}
      test_loc_1.p = location_positions[loc1_id];
      test_loc_2.p = location_positions[loc2_id];
      for (int vec1_id = 1; vec1_id < velocity_samples_tuples.size(); vec1_id++) {
        Vector<3> loc_1_velocity = std::get<0>(velocity_samples_tuples[vec1_id]);
        for (int vec2_id = 1; vec2_id < velocity_samples_tuples.size(); vec2_id++) {
          Vector<3> loc_2_velocity = std::get<0>(velocity_samples_tuples[vec2_id]);

          // Cost calculation.
          test_loc_1.v = loc_1_velocity;
          test_loc_2.v = loc_2_velocity;
          PointMassTrajectory3D tr(test_loc_1, test_loc_2, max_acc_per_axis, true);
          if (!tr.exists()) {
            std::cout << "Not-existing" << std::endl;
          }
          travel_costs[loc1_id][loc2_id][vec1_id][vec2_id] = tr.time();
          if (tr.time() <= 0 || tr.time() > 1000) {
            std::cout << "time inside new thing -> " << tr.time() << std::endl;
          }
        }
      }
    }
  }

  // Precalculate cost from [start_location](start_velocity) to other positions.
  // This will be used when `start_location` is not sampled.
  test_loc_1.p = location_positions[0];
  test_loc_1.v = std::get<0>(velocity_samples_tuples[0]);
  for (int loc2_id = 1; loc2_id < location_positions.size(); loc2_id++) {
    test_loc_2.p = location_positions[loc2_id];
    for (int vec2_id = 1; vec2_id < velocity_samples_tuples.size(); vec2_id++) {
      test_loc_2.v = std::get<0>(velocity_samples_tuples[vec2_id]);
      PointMassTrajectory3D tr(test_loc_1, test_loc_2, max_acc_per_axis, true);
      if (!tr.exists()) {
        std::cout << "Not-existing" << std::endl;
      }
      travel_costs[0][loc2_id][0][vec2_id] = tr.time();
      if (tr.time() <= 0 || tr.time() > 1000) {
        std::cout << "time inside new thing -> " << tr.time() << std::endl;
      }
    }
  }
  std::cout << "Finished." << std::endl;
};


void populate_precalculated_travel_costs_map(travel_cost_map &travel_costs,
                                             std::vector<Vector<3>> &location_positions,
                                             std::vector<Scalar> &velocity_norm_samples,
                                             std::vector<Scalar> &heading_angle_samples,
                                             norm_angle_to_velocity_vector_map &norm_angle_to_vector,
                                             Vector<3> &max_acc_per_axis,
                                             Vector<3> start_vel,
                                             Vector<3> end_vel) {


  std::tuple<Scalar, Scalar> start_vel_mag_ang = to_velocity_norm_and_angle(start_vel);
  QuadState test_loc_1;
  QuadState test_loc_2;
  test_loc_1.setZero();
  test_loc_2.setZero();
  for (int loc1_id = 0; loc1_id < location_positions.size(); loc1_id++) {
    for (int loc2_id = 0; loc2_id < location_positions.size(); loc2_id++) {
      if (loc1_id == loc2_id) { continue; }
      test_loc_1.p = location_positions[loc1_id];
      test_loc_2.p = location_positions[loc2_id];

      // For every combination of NORM and ANGLE
      for (Scalar norm1 : velocity_norm_samples) {
        for (Scalar norm2 : velocity_norm_samples) {
          for (Scalar angle1 : heading_angle_samples) {
            for (Scalar angle2 : heading_angle_samples) {
              Vector<3> loc_1_velocity = norm_angle_to_vector[norm1][angle1];
              Vector<3> loc_2_velocity = norm_angle_to_vector[norm2][angle2];

              // Location is START or END, the velocity is fixed (not sampled)
              if (loc1_id == 0) {
//                loc_1_velocity = start_vel;
//                norm1 = std::get<0>(start_vel_mag_ang);
//                angle1 = std::get<1>(start_vel_mag_ang);

                test_loc_1.v = start_vel;
                test_loc_2.v = loc_2_velocity;
                PointMassTrajectory3D tr(test_loc_1, test_loc_2, max_acc_per_axis, true);
                if (!tr.exists()) {
                  std::cout << "Not-existing" << std::endl;
                }
                travel_costs[loc1_id][loc2_id][std::get<0>(start_vel_mag_ang)][norm2][std::get<1>(start_vel_mag_ang)][angle2] = tr.time();
              }
              if (loc2_id == 0) {
//                loc_2_velocity = start_vel;
//                norm2 = std::get<0>(start_vel_mag_ang);
//                angle2 = std::get<1>(start_vel_mag_ang);

                test_loc_1.v = loc_1_velocity;
                test_loc_2.v = start_vel;
                PointMassTrajectory3D tr(test_loc_1, test_loc_2, max_acc_per_axis, true);
//                std::cout << "max_acc_per_axis: " << max_acc_per_axis << std::endl;
//                exit(1);
                if (!tr.exists()) {
                  std::cout << "Not-existing" << std::endl;
                }
                travel_costs[loc1_id][loc2_id][norm1][std::get<0>(start_vel_mag_ang)][angle1][std::get<1>(start_vel_mag_ang)] = tr.time();
              }

              test_loc_1.v = loc_1_velocity;
              test_loc_2.v = loc_2_velocity;
              PointMassTrajectory3D tr(test_loc_1, test_loc_2, max_acc_per_axis, true);
              if (!tr.exists()) {
                std::cout << "Not-existing" << std::endl;
              }
              //              if (loc1_id == 0 && loc2_id == 1) {
              //                std::cout << "from_p -> " << location_positions[loc1_id].transpose() << std::endl;
              //                std::cout << "to_p -> " << location_positions[loc2_id].transpose() << std::endl;
              //                std::cout << "from_v -> " << to_velocity_vector(norm1, angle1).transpose() << std::endl;
              //                std::cout << "to_v -> " << to_velocity_vector(norm2, angle2).transpose() << std::endl;
              //                std::cout << "Time between " <<loc1_id << " and " << loc2_id << " is " << tr.time() << std::endl;
              //              }
              travel_costs[loc1_id][loc2_id][norm1][norm2][angle1][angle2] = tr.time();
            }
          }
        }
      }
    }
  }

  //  auto time_taken = travel_costs[0][6][std::get<0>(start_vel_mag_ang)][velocity_norm_samples[1]][std::get<1>(start_vel_mag_ang)][heading_angle_samples[2]];
  //  std::cout << "From start to end -> " << time_taken << std::endl;
  //  exit(1);
}

std::vector<int> get_missing_values_in_range(std::vector<int> A, int start, int end) {
  std::sort(A.begin(), A.end());

  // Initialize the output vector
  std::vector<int> missing_nums;

  // Iterate over the range [start, end] and check for missing numbers
  for (int i = start; i <= end; i++) {
    // Use binary search to check if i is present in A
    if (!std::binary_search(A.begin(), A.end(), i)) {
      missing_nums.push_back(i);
    }
  }

  // Return the output vector
  return missing_nums;
}

std::vector<Scalar> get_mwp_trajectory_velocities(MultiWaypointTrajectory& trajectories) {
  std::vector<Scalar> velocities{};
  for (auto tr: trajectories) {
    velocities.push_back(tr.inp_from_v_norm);
  }
  velocities.push_back(trajectories[trajectories.size()-1].inp_to_v_norm);
  return velocities;
}

std::vector<Scalar> get_mwp_trajectory_yaw_angles(MultiWaypointTrajectory& trajectories) {
  std::vector<Scalar> yaw_angles{};
  for (auto tr: trajectories) {
    int ang1 = tr.inp_from_v_angle * 180 / M_PI;
    yaw_angles.push_back(ang1);
  }
  yaw_angles.push_back((int)(trajectories[trajectories.size()-1].inp_to_v_angle * 180 / M_PI));
  return yaw_angles;
}

Scalar get_mwp_trajectory_reward(std::vector<int>& scheduled_locations_idx, std::vector<Scalar>& rewards) {
  Scalar collected_reward = 0;
  for (int i = 1; i < scheduled_locations_idx.size(); i++) {
    collected_reward += rewards[scheduled_locations_idx[i]];
  }
  return collected_reward;
}

Scalar get_mwp_trajectory_cost(MultiWaypointTrajectory& trajectories) {
  Scalar cost = 0;
  for (auto& tr: trajectories) {
    cost += tr.time();
  }
  return cost;
}


} // namespace agi