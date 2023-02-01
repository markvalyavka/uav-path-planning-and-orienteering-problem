#pragma once

#include "gravity.hpp"
#include "pmm_trajectory3d.hpp"
#include "three_acc.hpp"
#include "timer.hpp"
#include "tuples_hash.hpp"
#include "velocity_search_graph.hpp"

#include "yaml-cpp/yaml.h"
#include "misc_helpers.hpp"

namespace agi {

class EnvConfig {
 public:
  EnvConfig(std::string cfg_file);

  void generate_samples_with_simple_sampling();
  void generate_precalculated_graph_of_costs();

  // FROM CFG YAML FILE
  Vector<3> start_velocity;
  Vector<3> start_position;
  Vector<3> end_velocity;
  Vector<3> end_position;
  int V;
  int H;
  Scalar t_max;
  Scalar max_velocity;
  Scalar max_acceleration;
  Vector<3> max_acc_per_axis;
  // Vector of 'heading angles'(positions) of all gates (including start and end).
  std::vector<Vector<3>> location_positions;
  std::vector<Scalar> rewards;

  // PRECALCULATED
  std::vector<Scalar> velocity_norm_samples;
  std::vector<Scalar> heading_angle_samples;
  std::vector<std::tuple<Vector<3>, Scalar, Scalar>> velocity_samples_tuples;
  norm_angle_to_velocity_vector_map norm_angle_to_vector;
  travel_cost_map precalculated_costs;
};

} // agi