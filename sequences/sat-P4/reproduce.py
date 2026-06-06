#!/usr/bin/env python3
"""Reproduce and verify sat(n, P_4): minimum edges in a P_4-saturated graph.

Generator: src/graph_saturation.py
  - sat_path(n, 4): exact MaxSAT (RC2) computation of sat(n, P_4).
  - sat_formula_p4(n): the proven Kaszonyi-Tuza (1986) closed form,
        a(n) = n/2 (n even), (n+3)/2 (n odd), for n >= 4.

Strategy (kept under ~60s):
  * FAST-RECOMPUTE sat(n, P_4) by MaxSAT for n = 4..9 (small instances).
  * Verify the ENTIRE published DATA (n = 4..23) against the proven
    closed form sat_formula_p4, which graph_saturation.py exposes and
    which the project validated against the MaxSAT solver for every n.

Run as:  python3 reproduce.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.graph_saturation import sat_path, sat_formula_p4  # noqa: E402

ID = "sat-P4"
OFFSET = 4
# Published DATA, a(n) for n = 4, 5, 6, ... (20 terms, n = 4..23).
DATA = [2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9, 8, 10, 9, 11, 10, 12, 11, 13]

# Recompute the small terms by MaxSAT (fast: n <= 9).  Larger terms are
# verified against the proven closed form only (MaxSAT on big n is slow).
FAST_MAX_N = 9


def main() -> None:
    n_terms = len(DATA)

    # (1) Verify every published term against the Kaszonyi-Tuza closed form
    #     that graph_saturation.py implements.
    for i, expected in enumerate(DATA):
        n = OFFSET + i
        formula = sat_formula_p4(n)
        assert formula == expected, (
            f"closed form sat_formula_p4({n}) = {formula} != DATA {expected}")

    # (2) Independently FAST-RECOMPUTE the small terms by MaxSAT and check
    #     they agree with both the closed form and DATA.
    n_sat = 0
    for i, expected in enumerate(DATA):
        n = OFFSET + i
        if n > FAST_MAX_N:
            break
        computed = sat_path(n, 4)
        assert computed == expected, (
            f"MaxSAT sat({n}, P_4) = {computed} != DATA {expected}")
        n_sat += 1

    print(f"OK {ID}: reproduced {n_terms} terms "
          f"({n_sat} via MaxSAT for n=4..{OFFSET + n_sat - 1}, "
          f"all {n_terms} via the proven closed form)")


if __name__ == "__main__":
    main()
