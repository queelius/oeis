#!/usr/bin/env python3
"""Reproduce and verify max-tail-bivariate.

The joint distribution B(n,t,c) = #{ endofunctions f:[n]->[n] : max tail
length (rho-height) T(f) = t AND f has exactly c cyclic points }, read by
rows as a full-square triangle (outer t = 0..n-1, inner c = 1..n, explicit
zeros), flattened over n = 1, 2, 3, ...

Generator: src/max_tail_bivariate.py
  - joint_tc_via_egf(n): exact counts via the proved bivariate EGF
        F_{<=t}(x,u) = 1/(1 - u e_t(x)),  B = n![x^n u^c](F_t - F_{t-1}).
  - brute_force_joint_tc(n): ground truth over all n^n maps.
  - cyclic_points_marginal(n): closed form for the sum-over-t marginal.
Marginals are cross-checked against src/max_tail_length.py
(distribution_triangle = OEIS A216242).

This reproduce.py (kept well under ~60s):
  * RECOMPUTES B via the EGF for n = 1..6 and checks it CELL-BY-CELL against
    brute force over all n^n maps.
  * VERIFIES both marginals exactly: sum_c B = A216242 row (from
    src/max_tail_length.distribution_triangle); sum_t B = A066324 row
    (cyclic_points_marginal).
  * BUILDS the flattened full-square DATA and asserts it equals the known
    DATA (n = 1..6, 91 terms).

Run as:  python3 reproduce.py
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.max_tail_bivariate import (  # noqa: E402
    brute_force_joint_tc,
    cyclic_points_marginal,
    joint_tc_rows,
    joint_tc_via_egf,
)
from src.max_tail_length import distribution_triangle  # noqa: E402

ID = "max-tail-bivariate"
N_MAX = 6  # brute force is feasible to n=6 (6^6 = 46656 maps)

# Known flattened full-square DATA B(n,t,c), n = 1..6 (outer t=0..n-1,
# inner c=1..n with explicit zeros).  See bivariate-README.md.
DATA = ("1,0,2,2,0,0,0,6,3,12,0,6,0,0,0,0,0,24,4,48,72,0,36,48,0,0,24,0,0,0,"
        "0,0,0,0,120,5,160,540,480,0,200,600,360,0,0,300,240,0,0,0,120,0,0,0,"
        "0,0,0,0,0,0,720,6,480,3240,5760,3600,0,1170,6000,7560,2880,0,0,3360,"
        "5040,2160,0,0,0,2520,1440,0,0,0,0,720,0,0,0,0,0")
KNOWN = [int(x) for x in DATA.split(",")]


def main() -> None:
    # (1) EGF == brute force, cell by cell; (2) both marginals exact.
    for n in range(1, N_MAX + 1):
        egf = joint_tc_via_egf(n)
        brute = brute_force_joint_tc(n)
        assert egf == brute, f"EGF != brute force at n={n}"

        # sum over c == A216242 row (max-tail distribution from max_tail_length)
        d_over_c = Counter()
        for (t, c), v in egf.items():
            d_over_c[t] += v
        a216242_row = distribution_triangle(n)[n - 1]
        assert [d_over_c[t] for t in range(n)] == a216242_row, \
            f"sum_c B != A216242 row at n={n}"

        # sum over t == A066324 row (cyclic-points marginal closed form)
        d_over_t = Counter()
        for (t, c), v in egf.items():
            d_over_t[c] += v
        assert [d_over_t[c] for c in range(1, n + 1)] == \
            cyclic_points_marginal(n), f"sum_t B != A066324 row at n={n}"

    # (3) Flattened full-square DATA matches the known DATA.
    flat = []
    for n in range(1, N_MAX + 1):
        for row in joint_tc_rows(n, full_square=True):
            flat.extend(row)
    assert flat == KNOWN, "flattened B(n,t,c) != known DATA"

    print(f"OK {ID}: reproduced {len(flat)} terms (B(n,t,c) full-square, "
          f"n=1..{N_MAX}); EGF == brute force cell-by-cell; marginals match "
          f"A216242 (max_tail_length) and A066324")


if __name__ == "__main__":
    main()
