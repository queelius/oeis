#!/usr/bin/env python3
"""Reproduce and verify crn-k2-by-reactions.

a(n) = number of chemical reaction networks on 2 species, up to permutation
of the species (S_2), having exactly n reactions, where every complex is
mono- or bimolecular (total stoichiometry <= 2), no complex is the zero
complex, and both species are used.  This is a complete, finite sequence:
the admissible directed-reaction ground set has R = 20 elements, so
a(n) = 0 for n > 20.

Generator: src/reaction_networks_census.py
  Two independent routes agree on every term:
    (A) orbit_count_burnside(2, 2, n): Burnside / cycle-index count of
        S_2-orbits of n-reaction subsets (all species used).  Instant.
    (B) census(k=2, max_stoich=2, r_max=None): direct S_2-canonical
        enumeration of all 2^20 reaction subsets.  ~40-50 s; OPT-IN here
        via env var CRN_FULL_CENSUS=1 to keep the default run fast.

Run as:        python3 reproduce.py
Full census:   CRN_FULL_CENSUS=1 python3 reproduce.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.reaction_networks_census import (  # noqa: E402
    census,
    orbit_count_burnside,
)

ID = "crn-k2-by-reactions"
OFFSET = 1
# Published DATA a(n) for n = 1..20 (then 0).
DATA = [8, 99, 570, 2445, 7752, 19440, 38760, 63090, 83980, 92504,
        83980, 63090, 38760, 19440, 7752, 2445, 570, 100, 10, 1]
GRAND_TOTAL = 524796  # sum over all n


def main() -> None:
    n_terms = len(DATA)

    # (A) PRIMARY: recompute every term via Burnside (instant, independent).
    burnside = [orbit_count_burnside(2, 2, r) for r in range(1, n_terms + 1)]
    assert burnside == DATA, f"Burnside {burnside} != DATA {DATA}"
    assert sum(burnside) == GRAND_TOTAL, \
        f"sum {sum(burnside)} != {GRAND_TOTAL}"

    routes = "Burnside cycle-index"
    if os.environ.get("CRN_FULL_CENSUS") == "1":
        # (B) CROSS-CHECK: full 2^20 direct S_2-canonical census (~40-50 s).
        workers = min(4, os.cpu_count() or 1)
        res = census(k=2, max_stoich=2, r_max=None, workers=workers)
        direct = res.r_total_sequence()
        assert direct == DATA, f"direct census {direct} != DATA {DATA}"
        assert res.total == GRAND_TOTAL, \
            f"census total {res.total} != {GRAND_TOTAL}"
        routes += " + full 2^20 direct census"

    print(f"OK {ID}: reproduced {n_terms} terms via {routes} "
          f"(sum = {GRAND_TOTAL})")


if __name__ == "__main__":
    main()
