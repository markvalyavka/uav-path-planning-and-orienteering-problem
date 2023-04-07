#include "env_config.hpp"

namespace agi {

EnvConfig::EnvConfig(std::string cfg_file) {
  YAML::Node config = YAML::LoadFile(cfg_file);

  config["start"]["velocity"] >> start_velocity;
  config["start"]["position"] >> start_position;
  config["end"]["velocity"] >> end_velocity;
  config["end"]["position"] >> end_position;

  V = config["number_of_velocity_samples"].as<int>();
  H = config["number_of_angle_samples"].as<int>();
  t_max = config["t_max"].as<Scalar>();
  max_velocity = config["max_velocity_norm"].as<Scalar>();
  Scalar max_acceleration = config["max_acceleration"].as<Scalar>();
  max_acc_per_axis = Vector<3>::Constant(max_acceleration);
  avg_reward_over_runs = config["avg_reward_over_runs"].as<int>();

  if (!parseArrayParam<Vector<3>>(config, "locations", location_positions)) {
    std::cerr << "YAML:PARSER -> Can't load param 'locations'" << std::endl;
  }
  location_positions.insert(location_positions.begin(), start_position);
  location_positions.push_back(end_position);

  if (!parseArrayParam<Scalar>(config, "rewards", rewards)) {
    std::cerr << "Can't load param 'rewards'" << std::endl;
  }
  rewards.insert(rewards.begin(), 0);
  rewards.push_back(0);

//  std::cout << "done." << std::endl;
}

void EnvConfig::generate_samples_with_simple_sampling() {
  velocity_norm_samples = get_velocity_norm_samples(V, max_velocity);
  heading_angle_samples = get_heading_angle_samples(H);

  // Samples combinations in vector form. [norm][angle] -> *vector [x, y, z]*
  populate_norm_angle_to_velocity_vector_map(norm_angle_to_vector,
                                             velocity_norm_samples,
                                             heading_angle_samples);

  velocity_samples_tuples = generate_total_sample_tuples(norm_angle_to_vector,
                                 velocity_norm_samples,
                                 heading_angle_samples);
}
void EnvConfig::generate_precalculated_graph_of_costs() {
  populate_precalculated_travel_costs_map(precalculated_costs,
                                          location_positions,
                                          velocity_norm_samples,
                                          heading_angle_samples,
                                          norm_angle_to_vector,
                                          max_acc_per_axis,
                                          start_velocity,
                                          end_velocity);

  imp_velocity_samples_tuples = velocity_samples_tuples;
  imp_velocity_samples_tuples.insert(imp_velocity_samples_tuples.begin(), {start_velocity,
                                                                   std::get<0>(to_velocity_norm_and_angle(start_velocity)),
                                                                   std::get<1>(to_velocity_norm_and_angle(start_velocity))});
  imp_populate_precalculated_travel_costs_map(imp_precalculated_costs,
                                          location_positions,
                                          velocity_samples_tuples,
                                          max_acc_per_axis);
};

} // namespace agi