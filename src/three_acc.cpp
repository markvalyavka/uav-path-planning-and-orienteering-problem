#include "three_acc.hpp"

namespace agi {

void three_acc(const Scalar ps, const Scalar vs, const Scalar pe,
               const Scalar ve, const Scalar a1, const Scalar a2,
               const Scalar scale, const Scalar tswitch) {
  std::cout << "a1 " << a1 << " a2 " << a2 << " scale " << scale << std::endl;
  const Scalar pow_ve2 = ve * ve;
  const Scalar pow_vs2 = vs * vs;
  const Scalar pow_a1_2 = a1 * a1;
  const Scalar pow_a2_2 = a2 * a2;
  //   const Scalar pow_a3_2 = a3 * a3;
  const Scalar pow_tswitch_2 = tswitch * tswitch;

  Scalar sqrt_part =
    sqrt((a1 - a2) *
         (a1 * (2 * a2 * (-pe + ps) * scale +
                pow_a2_2 * (-1 + scale) * scale * pow_tswitch_2 + pow_ve2) -
          a2 * (pow_ve2 - scale * pow_ve2 + scale * pow_vs2 +
                2 * a2 * (-1 + scale) * scale * (pe - ps - tswitch * vs))));
  Scalar t1_1 = (pow_a2_2 * tswitch + a1 * a2 * (-1 + scale) * tswitch -
                 pow_a2_2 * scale * tswitch - a1 * vs + a2 * vs + sqrt_part) /
                ((a1 - a2) * (a1 + a2 * (-1 + scale)));
  Scalar t2_1 =
    (pow_a1_2 * tswitch - a2 * vs + a1 * (-(a2 * tswitch) + vs) - sqrt_part) /
    ((a1 - a2) * (a1 + a2 * (-1 + scale)));
  Scalar t3_1 = -((a1 * a2 * scale * tswitch - a1 * ve + a2 * ve -
                   a2 * scale * ve + a2 * scale * vs + sqrt_part) /
                  (a2 * (a1 + a2 * (-1 + scale)) * scale));
  Scalar t1_2 = (pow_a2_2 * tswitch + a1 * a2 * (-1 + scale) * tswitch -
                 pow_a2_2 * scale * tswitch - a1 * vs + a2 * vs - sqrt_part) /
                ((a1 - a2) * (a1 + a2 * (-1 + scale)));
  Scalar t2_2 =
    (pow_a1_2 * tswitch - a2 * vs + a1 * (-(a2 * tswitch) + vs) + sqrt_part) /
    ((a1 - a2) * (a1 + a2 * (-1 + scale)));
  Scalar t3_2 =
    (-(a2 * ve) + a2 * scale * ve + a1 * (-(a2 * scale * tswitch) + ve) -
     a2 * scale * vs + sqrt_part) /
    (a2 * (a1 + a2 * (-1 + scale)) * scale);

  //   Scalar sqrt_part =
  //     sqrt((a1 - a2) *
  //          (a1 * (pow_a3_2 * pow_tswitch_2 -
  //                 a3 * (2 * pe - 2 * ps + a2 * pow_tswitch_2) + pow_ve2) -
  //           a2 * (pow_ve2 + 2 * a3 * (-pe + ps + tswitch * vs)) +
  //           a3 * (pow_ve2 - pow_vs2 + 2 * a3 * (-pe + ps + tswitch * vs))));
  //   Scalar t1_1 = (-(a1 * a2 * tswitch) + pow_a2_2 * tswitch + a1 * a3 *
  //   tswitch -
  //                  a2 * a3 * tswitch - a1 * vs + a2 * vs + sqrt_part) /
  //                 ((a1 - a2) * (a1 - a2 + a3));
  //   Scalar t2_1 = -((-(pow_a1_2 * tswitch) + a1 * a2 * tswitch - a1 * vs +
  //                    a2 * vs + sqrt_part) /
  //                   ((a1 - a2) * (a1 - a2 + a3)));
  //   Scalar t3_1 =
  //     -((a1 * a3 * tswitch - a1 * ve + a2 * ve - a3 * ve + a3 * vs +
  //     sqrt_part) /
  //       (a3 * (a1 - a2 + a3)));

  //   Scalar t1_2 = -((a1 * a2 * tswitch - pow_a2_2 * tswitch - a1 * a3 *
  //   tswitch +
  //                    a2 * a3 * tswitch + a1 * vs - a2 * vs + sqrt_part) /
  //                   ((a1 - a2) * (a1 - a2 + a3)));
  //   Scalar t2_2 =
  //     (pow_a1_2 * tswitch - a1 * a2 * tswitch + a1 * vs - a2 * vs +
  //     sqrt_part) /
  //     ((a1 - a2) * (a1 - a2 + a3));
  //   Scalar t3_2 =
  //     (-(a1 * a3 * tswitch) + a1 * ve - a2 * ve + a3 * ve - a3 * vs +
  //     sqrt_part) / (a3 * (a1 - a2 + a3));

  std::cout << "t1_1 " << t1_1 << " t2_1 " << t2_1 << " t3_1 " << t3_1
            << std::endl;
  std::cout << "t1_2 " << t1_2 << " t2_2 " << t2_2 << " t3_2 " << t3_2
            << std::endl;
}

}  // namespace agi