#pragma once

#include <cmath>
#include <iostream>

#include "agilib/math/types.hpp"

namespace agi {

void three_acc(const Scalar ps, const Scalar vs, const Scalar pe,
               const Scalar ve, const Scalar a1, const Scalar a2,
               const Scalar scale, const Scalar tswitch);

}  // namespace agi