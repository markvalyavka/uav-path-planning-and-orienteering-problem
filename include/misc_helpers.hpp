#pragma once


#include <unordered_map>
#include "gravity.hpp"
#include "pmm_trajectory3d.hpp"
#include "three_acc.hpp"
#include "timer.hpp"
#include "tuples_hash.hpp"
#include "velocity_search_graph.hpp"
#include "yaml-cpp/yaml.h"




namespace agi {

typedef std::unordered_map<int, std::unordered_map<int,
                                                   std::unordered_map<Scalar,std::unordered_map<Scalar,
                                                                                                 std::unordered_map<Scalar, std::unordered_map<Scalar, Scalar>>>>>>
  travel_cost_map;

typedef std::unordered_map<Scalar, std::unordered_map<Scalar, Vector<3>>>
  norm_angle_to_velocity_vector_map;

std::vector<Scalar> get_velocity_norm_samples(int samples_num, Scalar max_velocity);
std::vector<Scalar> get_heading_angle_samples(int samples_num);

/*
 * Returns velocity vector represented as <X, Y, Z> coordinates given
 * `velocity_norm` and `angle`.
 */
Vector<3> to_velocity_vector(Scalar velocity_norm, Scalar angle);

/*
 * Returns velocity vector represented as magnitude(norm) and direction(angle)
 * given <X, Y, Z> coordinates.
 */
std::tuple<Scalar, Scalar> to_velocity_norm_and_angle(Vector<3> velocity_vector);


std::vector<std::tuple<Vector<3>, Scalar, Scalar>> generate_total_sample_tuples(
  norm_angle_to_velocity_vector_map &norm_angle_to_vector,
  std::vector<Scalar> &velocity_norm_samples,
  std::vector<Scalar> &heading_angle_samples);


void populate_norm_angle_to_velocity_vector_map(norm_angle_to_velocity_vector_map &map,
                                                std::vector<Scalar> &vel_norm_samples,
                                                std::vector<Scalar> &angle_samples);

void populate_precalculated_travel_costs_map(travel_cost_map &travel_costs,
                                             std::vector<Vector<3>> &location_positions,
                                             std::vector<Scalar> &velocity_norm_samples,
                                             std::vector<Scalar> &heading_angle_samples,
                                             norm_angle_to_velocity_vector_map &norm_angle_to_vector,
                                             Vector<3> &max_acc_per_axis,
                                             Vector<3> start_vel,
                                             Vector<3> end_vel);


} // namespace agi