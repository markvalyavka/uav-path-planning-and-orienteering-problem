#pragma once
#include <vector>

#include "pmm_trajectory3d.hpp"
#include "yaml-cpp/yaml.h"
namespace agi {

void operator>>(const YAML::Node& node, Scalar& value);
void operator>>(const YAML::Node& node, agi::Vector<3>& value);

template<typename DataType>
bool parseArrayParam(const YAML::Node& config, std::string param,
                     std::vector<DataType>& array) {
  if (config[param.c_str()]) {
    YAML::Node vector = config[param.c_str()];
    for (YAML::const_iterator ti = vector.begin(); ti != vector.end(); ++ti) {
      const YAML::Node& datayaml = *ti;
      DataType from_yaml;
      datayaml >> from_yaml;
      array.push_back(from_yaml);
    }
  } else {
    return false;
  }
  return true;
}

typedef std::vector<PointMassTrajectory3D> MultiWaypointTrajectory;

class VelocitySearchGraph {
 public:
  VelocitySearchGraph(const Scalar max_acc, const Scalar max_yaw_pitch_ang,
                      const Scalar precision_yaw_pitch_ang,
                      const Scalar yaw_pitch_cone_angle_boundary,
                      const Scalar min_velocity_size_boundary,
                      const Scalar min_velocity_size,
                      const Scalar max_velocity_size,
                      const Scalar precision_velocity_size);


  MultiWaypointTrajectory find_velocities_in_positions(
    const std::vector<Vector<3>>& gates_waypoints,
    const Vector<3>& start_velocity, const Vector<3>& end_velocity,
    const std::vector<Scalar>& gates_yaw_deg,
    const std::vector<Scalar>& gates_pitch_deg,
    const std::vector<Scalar>& gates_vel_norms, const bool end_free,
    const bool use_gd);
  void setYawPitchConeSidesAndPecisiton(const Scalar max_yaw_pitch_ang,
                                        const Scalar precision_yaw_pitch_ang);

  static void saveTrajectoryEquitemporal(
    const MultiWaypointTrajectory& trajectories, std::string filename,
    const Scalar dt_desired = 0.01);
  static void saveTrajectoryEquidistant(
    const MultiWaypointTrajectory& trajectories, std::string filename,
    const Scalar desired_ds = 0.250);
  static std::vector<QuadState> getTrajectoryEquidistantStates(
    const MultiWaypointTrajectory& trajectories,
    const Scalar desired_ds = 0.250);

  static std::vector<Vector<3>> sampleTrajectory(
    const MultiWaypointTrajectory& trajectory, const Scalar ds_desired);
  static void saveSamplesToFile(std::string filename,
                                std::vector<QuadState> samples);
  static std::string stateToStringHeader(const QuadState& s);
  static std::string stateToString(const QuadState& s);

 private:
  const Scalar max_acc_;
  Scalar max_yaw_pitch_ang_;
  Scalar precision_yaw_pitch_ang_;
  const Scalar max_velocity_size_;
  const Scalar min_velocity_size_;
  const Scalar min_velocity_size_boundary_;
  const Scalar yaw_pitch_cone_angle_boundary_;
  const Scalar precision_velocity_size_;
};

}  // namespace agi