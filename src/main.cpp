
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <csignal>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_set>

#include "agilib/math/gravity.hpp"
#include "agilib/utils/timer.hpp"
#include "pmm_trajectory3d.hpp"
#include "three_acc.hpp"
#include "tuples_hash.hpp"
#include "velocity_search_graph.hpp"
#include "yaml-cpp/yaml.h"

using namespace agi;

void signal_callback(int sig) {
  std::cout << "Signal " << sig << " received" << std::endl;
  exit(sig);
}

struct node {
  QuadState state;
  Scalar g = 0;  // current cost
  Scalar h = 0;  // heuristic
  int sample_index;
  node* parent;
  Scalar f() { return g + h; };
};

typedef std::tuple<short, short, short, short, short, short, short, short,
                   short>
  dist_state;

dist_state make_diststate(QuadState state) {
  return std::make_tuple((short)(state.p(0) * 5.0), (short)(state.p(1) * 5.0),
                         (short)(state.p(2) * 5.0), (short)(state.v(0) * 2.0),
                         (short)(state.v(1) * 2.0), (short)(state.v(2) * 2.0),
                         (short)(state.a(0) / 3.0), (short)(state.a(1) / 3.0),
                         (short)(state.a(2) / 3.0));
}

node* test_astar(MultiWaypointTrajectory tr) {
  auto cmp_node = [](node* l, node* r) { return l->f() > r->f(); };
  std::priority_queue<node*, std::vector<node*>, decltype(cmp_node)> open_set(
    cmp_node);  // priority queue for a*
  std::unordered_map<dist_state, node*> open_set_graph;

  std::unordered_set<dist_state> close_set;
  std::vector<QuadState> samples_dense =
    VelocitySearchGraph::getTrajectoryEquidistantStates(tr, 0.8);
  Scalar time_to_end = samples_dense.back().t;
  node* root = new node();
  root->state = samples_dense[0];
  root->state.a = Vector<3>(0, 0, 0);
  root->g = 0;
  root->h = time_to_end;
  root->sample_index = 0;
  root->parent = NULL;
  open_set.push(root);

  auto root_dist = make_diststate(root->state);
  open_set_graph[root_dist] = root;
  close_set.insert(root_dist);

  Scalar ds = 0.4;
  Scalar max_thrust = 28;
  int iter = 0;
  Scalar furthest = 0;
  node* solution = NULL;
  while (open_set.size() > 0) {
    iter++;
    node* cbest = open_set.top();
    // std::cout << "cbest f " << cbest->f() << " time " << cbest->state.t
    //           << std::endl;
    if (cbest->g > furthest) {
      furthest = cbest->g;
      std::cout << " new best f " << cbest->f() << " time " << cbest->state.t
                << std::endl;
    }
    Scalar dist = (cbest->state.p - samples_dense.back().p).norm();
    // std::cout << "best dist " << dist << std::endl;
    if (dist < 0.4) {
      std::cout << "reached final state" << std::endl;
      solution = cbest;
      break;
    }
    open_set.pop();
    dist_state cbest_dist = make_diststate(cbest->state);
    close_set.insert(cbest_dist);

    Scalar t = ds / cbest->state.v.norm();
    // Scalar t = 0.1;
    Scalar t_a = sqrt(2 * ds / max_thrust);
    if (!std::isfinite(t)) {
      t = t_a;
    }

    // std::cout << "sim for " << t << std::endl;
    Scalar t_pow2 = t * t;
    Scalar t_pow3 = t_pow2 * t;
    QuadState cbest_sample = samples_dense[cbest->sample_index];
    QuadState cbest_sample_next = samples_dense[std::min(
      cbest->sample_index + 1, (int)samples_dense.size() - 1)];

    for (auto jx : {300.0, 150.0, 0.0, -150.0, -300.0}) {
      for (auto jy : {300.0, 150.0, 0.0, -150.0, -300.0}) {
        for (auto jz : {300.0, 150.0, 0.0, -150.0, -300.0}) {
          Vector<3> j(jx, jy, jz);
          Vector<3> new_a = cbest->state.a + j * t;
          Vector<3> thrust = new_a - GVEC;
          if (thrust.norm() > max_thrust) {
            // std::cout << "above" << std::endl;
            continue;
          }


          Vector<3> new_p = cbest->state.p + cbest->state.v * t +
                            0.5 * cbest->state.a * t_pow2 +
                            (1.0 / 6.0) * j * t_pow3;
          Vector<3> new_v =
            cbest->state.v + cbest->state.a * t + 0.5 * j * t_pow2;
          Scalar dist_cur = (new_p - cbest_sample.p).norm();
          Scalar dist_next = (new_p - cbest_sample_next.p).norm();
          Scalar dist_cur_v = (new_v - cbest_sample.v).norm();
          Scalar dist_next_v = (new_v - cbest_sample_next.v).norm();

          if (dist_cur + dist_next > 1.5) {
            // std::cout << "too far" << std::endl;
            continue;
          }
          if (dist_cur_v + dist_next_v > 8.0) {
            // std::cout << "too far" << std::endl;
            continue;
          }
          // std::cout << "ceck if the expanded one is in closed list!!!!"
          //           << std::endl;
          // exit(1);

          Scalar new_g = cbest->g + t;
          QuadState new_state;
          new_state.p = new_p;
          new_state.v = new_v;
          new_state.a = new_a;
          new_state.j = j;
          new_state.t = cbest->g + t;
          dist_state new_state_dist = make_diststate(new_state);

          node* newnode;
          auto search_closed = close_set.find(new_state_dist);

          bool add_to_open = false;
          if (search_closed == close_set.end()) {  // not in closed

            // newer is better
            auto search_open = open_set_graph.find(new_state_dist);
            if (search_open == open_set_graph.end()) {
              // not in open set
              add_to_open = true;
            } else {
              if (new_g > search_open->second->g) {
                // is in open set but already better == do not consider further
                continue;
              } else {
                newnode = search_open->second;
              }
            }
          }


          if (add_to_open) {
            // has not been there
            newnode = new node();
          }

          newnode->state = new_state;
          newnode->g = new_g;
          newnode->h = (time_to_end - cbest_sample.t -
                        dist_cur / (dist_cur + dist_next) *
                          (cbest_sample_next.t - cbest_sample.t));

          // move sample index if appropriate
          newnode->sample_index = dist_cur > dist_next ? cbest->sample_index + 1
                                                       : cbest->sample_index;
          newnode->parent = cbest;
          // std::cout << "add node " << newnode->f() << " p "
          //           << new_state.p.transpose() << " v "
          //           << new_state.v.transpose() << " a "
          //           << new_state.a.transpose() << " j "
          //           << new_state.j.transpose() << std::endl;
          open_set.push(newnode);
          open_set_graph[new_state_dist] = newnode;
        }
      }
    }
    // if (iter >= 10) exit(1);
    // std::cout << "loop" << std::endl;
  }
  return solution;
}

int test_pmm(int argc, char** argv) {
  // register singal for killingcm
  std::signal(SIGINT, signal_callback);
  std::cout << "Testing PMM trajectories " << std::endl;

  std::string config_file = "config.yaml";
  YAML::Node config = YAML::LoadFile(config_file);

  QuadState from;
  from.setZero();
  from.p = Vector<3>(0, 0, 1);
  QuadState to;
  to.setZero();
  to.p = Vector<3>(5, 5, 1);
  PointMassTrajectory3D test(from, to, 25);
  std::cout << test << std::endl;

  // -+ value of initial samples
  const Scalar max_yaw_pitch_ang = 15.0;
  // distance between initial samples
  const Scalar precision_yaw_pitch_ang = 15.0;
  // maximal values of the samples
  const Scalar yaw_pitch_cone_angle_boundary = 60.0;

  // velocity norm samples
  const Scalar min_velocity_norm = 1.0;
  const Scalar min_velocity_norm_boundary = 1.0;
  const Scalar max_velocity_norm = 17.0;
  const Scalar precision_velocity_norm = 8.0;

  const Scalar max_acc_norm = 32.94;

  VelocitySearchGraph vel_search_graph(
    max_acc_norm, max_yaw_pitch_ang, precision_yaw_pitch_ang,
    yaw_pitch_cone_angle_boundary, min_velocity_norm_boundary,
    min_velocity_norm, max_velocity_norm, precision_velocity_norm);

  Vector<3> start_velocity;
  Vector<3> end_velocity;
  Vector<3> start_position;
  Vector<3> end_position;
  config["start"]["velocity"] >> start_velocity;
  config["end"]["velocity"] >> start_velocity;
  config["start"]["position"] >> start_position;
  config["end"]["position"] >> end_position;
  const bool end_free = true;
  std::vector<Vector<3>> gates_waypoints;
  std::vector<Scalar> gates_yaw_deg;

  if (!parseArrayParam<Vector<3>>(config, "gates", gates_waypoints))
    std::cerr << "can not load param gates" << std::endl;

  if (!parseArrayParam<Scalar>(config, "gates_orientations", gates_yaw_deg))
    std::cerr << "can not load param gates_orientations" << std::endl;

  gates_waypoints.insert(gates_waypoints.begin(), start_position);
  gates_waypoints.push_back(end_position);

  std::cout << "gates_waypoints.size() " << gates_waypoints.size() << std::endl;
  std::cout << "gates_yaw_deg.size() " << gates_yaw_deg.size() << std::endl;


  gates_waypoints.resize(3);
  gates_yaw_deg.resize(3);
  Scalar sum_times = 0;
  MultiWaypointTrajectory tr = vel_search_graph.find_velocities_in_positions(
    gates_waypoints, start_velocity, end_velocity, gates_yaw_deg, end_free,
    false);
  std::cout << "output tr size " << tr.size() << std::endl;
  for (size_t i = 0; i < tr.size(); i++) {
    if (i == 0) {
      std::cout << i << " vel " << tr[i].get_start_state().v.transpose()
                << " acc " << tr[i].get_start_state().a.transpose()
                << " thrust norm " << (tr[i].get_start_state().a - GVEC).norm()
                << std::endl;
    }
    std::cout << i + 1 << " vel " << tr[i].get_end_state().v.transpose()
              << " acc " << tr[i].get_end_state().a.transpose()
              << " thrust norm " << (tr[i].get_end_state().a - GVEC).norm()
              << std::endl;

    // std::cout << tr[i] << std::endl;
    sum_times += tr[i].time();
    if (tr[i].time() - tr[i].time_min() > PRECISION_PMM_VALUES) {
      std::cout << "bad time!!!!!!" << std::endl;
    }
  }


  PointMassTrajectory3D pmm3d = tr[0];
  std::cout << "pmm3d" << std::endl;
  std::cout << pmm3d << std::endl;
  QuadState from_pmm3d = pmm3d.get_start_state();
  QuadState to_pmm3d = pmm3d.get_end_state();

  Timer t_new("new pmm thrust");
  t_new.tic();
  for (size_t i = 0; i < 1000; i++) {
    PointMassTrajectory3D tr(from_pmm3d, to_pmm3d, max_acc_norm, 10);
  }
  t_new.toc();
  t_new.print();

  PointMassTrajectory3D tstpmm3d(from_pmm3d, to_pmm3d, max_acc_norm, 10);
  std::cout << "tstpmm3d " << tstpmm3d << std::endl;
  /*
  std::cout << "pmm start acc" << pmm3d.get_start_state().a.norm() << std::endl;
  std::cout << "pmm end acc" << pmm3d.get_end_state().a.norm() << std::endl;
  int iter = 0;
  while (fabs((pmm3d.get_end_state().a - GVEC).norm() - max_acc_norm) > 0.1) {
    iter++;
    Vector<3> acc_vec = pmm3d.get_start_state().a.cwiseAbs();
    Scalar acc_vec_norm = acc_vec.norm();
    // Vector<3> max_acc_new =
    //   acc_vec * sqrt((26.29 * 26.29) / (acc_vec_norm * acc_vec_norm));
    Scalar a_dis = acc_vec_norm * acc_vec_norm;
    Scalar b_dis = 2 * acc_vec(2) * G;
    Scalar c_dis = G * G - (max_acc_norm * max_acc_norm);
    Scalar k1 =
      (-b_dis + sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    Scalar k2 =
      (-b_dis - sqrt(b_dis * b_dis - 4 * a_dis * c_dis)) / (2 * a_dis);
    Vector<3> thrust_acc_new_k1 = k1 * acc_vec - GVEC;
    Vector<3> thrust_acc_new_k2 = k2 * acc_vec - GVEC;
    Vector<3> max_acc_new1 = thrust_acc_new_k1 + GVEC;
    Vector<3> max_acc_new2 = -thrust_acc_new_k1 + GVEC;
    std::cout << "acc_vec " << acc_vec.transpose() << std::endl;
    std::cout << "max_acc_new1 " << max_acc_new1.transpose() << std::endl;
    std::cout << "max_acc_new2 " << max_acc_new2.transpose() << std::endl;
    std::cout << "k1 " << k1 << " k1 * acc_vec " << (k1 * acc_vec).transpose()
              << " thrust_acc_new_k1 " << thrust_acc_new_k1.transpose()
              << " norm " << thrust_acc_new_k1.norm() << std::endl;
    std::cout << "k2 " << k2 << " k2 * acc_vec " << (k2 * acc_vec).transpose()
              << " thrust_acc_new_k2 " << thrust_acc_new_k2.transpose()
              << " norm " << thrust_acc_new_k2.norm() << std::endl;

    pmm3d =
      PointMassTrajectory3D(from_pmm3d, to_pmm3d, max_acc_new1, max_acc_new2);
    std::cout << "pmm3d iter " << iter << std::endl;
    std::cout << pmm3d << std::endl;
    std::cout << "pmm start acc" << pmm3d.get_start_state().a.norm()
              << " thrust norm " << (pmm3d.get_start_state().a - GVEC).norm()
              << std::endl;
    std::cout << "pmm end acc" << pmm3d.get_end_state().a.norm()
              << " thrust norm " << (pmm3d.get_end_state().a - GVEC).norm()
              << std::endl;
  }

*/
  // PointMassTrajectory3D pmm_gd(from_pmm3d, to_pmm3d, 26.29);
  // std::cout << "pmm_gd " << pmm_gd << std::endl;

  Timer t_gd("pmm gd");
  t_gd.tic();
  for (size_t i = 0; i < 1000; i++) {
    PointMassTrajectory3D tr(from_pmm3d, to_pmm3d, max_acc_norm);
  }
  t_gd.toc();
  t_gd.print();

  PointMassTrajectory3D pmm_gd2(from_pmm3d, to_pmm3d, max_acc_norm);
  std::cout << "pmm_gd2 " << pmm_gd2 << std::endl;
  // std::cout << tr[0] << std::endl;

  // std::cout << "case 1" << std::endl;
  // three_acc(tr[0].y_.p_(0), tr[0].y_.v_(0), tr[0].y_.p_(3), tr[0].y_.v_(2),
  //           tr[0].y_.a_(0), tr[0].y_.a_(1), 1.5, tr[0].z_.t_(0));
  // std::cout << "case 2" << std::endl;
  // three_acc(tr[0].y_.p_(0), tr[0].y_.v_(0), tr[0].y_.p_(3), tr[0].y_.v_(2),
  //           tr[0].y_.a_(1), tr[0].y_.a_(0), 1.5, tr[0].z_.t_(0));
  // std::cout << "case 3" << std::endl;
  // three_acc(tr[0].y_.p_(0), tr[0].y_.v_(0), tr[0].y_.p_(3), tr[0].y_.v_(2),
  //           tr[0].y_.a_(0), tr[0].y_.a_(1), -1.5, tr[0].z_.t_(0));
  // std::cout << "case 4" << std::endl;
  // three_acc(tr[0].y_.p_(0), tr[0].y_.v_(0), tr[0].y_.p_(3), tr[0].y_.v_(2),
  //           tr[0].y_.a_(1), tr[0].y_.a_(0), -1.5, tr[0].z_.t_(0));

  // node* solution = test_astar(tr);

  // std::vector<node*> traj;
  // node* cur = solution;
  // while (cur != NULL) {
  //   traj.push_back(cur);
  //   cur = cur->parent;
  // }
  // if (traj.size() > 0) {
  //   std::reverse(traj.begin(), traj.end());

  //   std::cout << "trajectory:" << std::endl;
  //   for (size_t i = 0; i < traj.size(); i++) {
  //     std::cout << traj[i]->state.p.transpose() << std::endl;
  //   }

  //   std::cout << "total time sampling " << traj.back()->state.t << std::endl;
  // }

  std::cout << "total time pmm " << sum_times << std::endl;
  std::vector<QuadState> samples =
    VelocitySearchGraph::getTrajectoryEquidistantStates(tr, 0.8);

  // A*x = b
  // x coefficients
  // A polynomial coefficients
  // b = [p0,p1,v0,v1,a0,a1]
  std::vector<std::vector<Vector<6>>> polynomials;
  std::cout << "samples.size() " << samples.size() << std::endl;
  polynomials.resize(samples.size() - 1);
  for (size_t i = 1; i < samples.size(); i++) {
    QuadState& from = samples[i - 1];
    QuadState& to = samples[i];
    Scalar dt = to.t - from.t;
    // x = [a5,a4,a3,a2,a1,a0]
    for (size_t axi = 0; axi < 3; axi++) {
      Vector<8> b;
      b << from.p(axi), to.p(axi), from.v(axi), to.v(axi), from.a(axi),
        to.a(axi), 0, 0;
      Matrix<8, 6> A;
      Vector<6> tau = Vector<6>::Ones();
      for (int i = 1; i < 6; i++) {
        tau(i) = tau(i - 1) * dt;
      }
      // std::cout << "tau " << tau.transpose() << std::endl;
      A.row(0) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;                    // p0
      A.row(1) << tau(5), tau(4), tau(3), tau(2), tau(1), tau(0);  // p1
      A.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;                    // v0
      A.row(3) << 5.0 * tau(4), 4.0 * tau(3), 3.0 * tau(2), 2.0 * tau(1),
        1.0 * tau(0), 0.0;                       // v1
      A.row(4) << 0.0, 0.0, 0.0, 2.0, 0.0, 0.0;  // a0
      A.row(5) << 20.0 * tau(3), 12.0 * tau(2), 6.0 * tau(1), 2.0 * tau(0), 0.0,
        0.0;                                     // a1
      A.row(6) << 0.0, 0.0, 6.0, 0.0, 0.0, 0.0;  // j0
      A.row(7) << 60.0 * tau(2), 24.0 * tau(1), 6.0 * tau(0), 0, 0.0,
        0.0;  // j1
      // A.row(6) << 0.0, 24.0, 0.0, 0.0, 0.0, 0.0;  // j0
      // A.row(7) << 120.0 * tau(1), 24.0 * tau(0), 0.0, 0.0, 0.0,
      //   0.0;  // j1

      // std::cout << "solve" << std::endl;
      Vector<6> p = A.colPivHouseholderQr().solve(b);
      polynomials[i - 1].push_back(p);
      // std::cout << "solved" << std::endl;
      // std::cout << "p" << i << "[" << axi << "] " << p.transpose() <<
      // std::endl;
    }
  }

  std::ofstream myfile;
  myfile.open("polynomials.csv");
  if (myfile.is_open()) {
    if (samples.size() > 0 && samples[0].size() > 0) {
      myfile << "i,axi,tfrom,tto,a5,a4,a3,a2,a1,a0" << std::endl;
    }

    for (int var1 = 0; var1 < polynomials.size(); ++var1) {
      const Scalar tfrom = samples[var1].t;
      const Scalar tto = samples[var1 + 1].t;
      for (size_t axi = 0; axi < 3; axi++) {
        const Vector<6>& p = polynomials[var1][axi];
        myfile << var1 << "," << axi << "," << tfrom << "," << tto << ","
               << p(0) << "," << p(1) << "," << p(2) << "," << p(3) << ","
               << p(4) << "," << p(5);
        myfile << std::endl;
      }
    }

    myfile.close();
  }

  std::cout << "out" << std::endl;

  VelocitySearchGraph::saveTrajectoryEquitemporal(tr, "samples_pmm.csv");

  std::cout << "saved equitemporal" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(tr, "samples_equidistant.csv");
  std::cout << "saved equidistant" << std::endl;
  VelocitySearchGraph::saveTrajectoryEquidistant(
    tr, "samples_equidistant_08.csv", 0.8);
  std::cout << "saved equidistant 0.8" << std::endl;

  return 0;
}

int main(int argc, char** argv) {
  test_pmm(argc, argv);
  return 0;
}
