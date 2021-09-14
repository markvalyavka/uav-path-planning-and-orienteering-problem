#include "pmm_trajectory3d.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace agi {

PointMassTrajectory3D::PointMassTrajectory3D(const QuadState &from,
                                             const QuadState &to,
                                             Scalar max_acc) {
  x_ =
    PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc, max_acc, 0);
  y_ =
    PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc, max_acc, 1);
  z_ =
    PMMTrajectory(from.p(0), from.v(0), to.p(0), to.v(0), max_acc, max_acc, 2);
}

std::ostream &operator<<(std::ostream &o, const PointMassTrajectory3D &t) {
  o << "pmm3d: t:" << f.time() << ";exists:" << f.exists();
  o << "\n\tx: " << t.x;
  o << "\n\ty: " << t.y;
  o << "\n\tz: " << t.z;
  return o;
}

}  // namespace agi