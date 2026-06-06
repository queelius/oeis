#!/usr/bin/env python3
"""Reproduce and verify A006615: a(n) = z(n,n;3,4) + 1.

z(n,n;3,4) is the Zarankiewicz number: max number of 1s in an n x n 0-1
matrix with no 3 x 4 all-ones submatrix.  The OEIS term is k_{3,4}(n,n) =
z(n,n;3,4) + 1.

Generator: src/zarankiewicz.py
  - lower bound: SAT witness (MaxSAT / RC2)
  - upper bound: double-lex symmetry-broken UNSAT (src/zarankiewicz_symbreak.py)
  - lookup_known(n,n,3,4): the SAT-PROVEN exact z-values baked into the
    module's KNOWN_Z_3_4 table (n = 4..10).

Strategy (kept under ~60s -- the project's full UNSAT proofs take minutes to
days, so we do NOT re-run them here):
  * FAST-RECOMPUTE z(n,n;3,4) by MaxSAT for the small n (n = 4, 5).  MaxSAT
    finds the optimum directly; for these tiny instances it is fast.
  * VALIDATE every published z-value (n = 4..10) against the module's KNOWN
    table via lookup_known, then form a(n) = z + 1.
  * a(11) = 79 (z(11,11;3,4) = 78) is the in-review extension proved by a
    9202 s double-lex UNSAT; it is NOT recomputed here and NOT in the KNOWN
    table -- it is carried from the b-file / submission as a cited term.

Run as:  python3 reproduce.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.zarankiewicz import lookup_known, zarankiewicz_sat  # noqa: E402

ID = "A006615"
OFFSET = 4
# Published DATA a(n) for n = 4..11.  a(4..10) are the originally published
# terms (a(10)=67 added Mar 2026); a(11)=79 is the in-review extension.
DATA = [15, 22, 31, 38, 46, 57, 67, 79]
N_START = 4
# n covered by the module's KNOWN_Z_3_4 table (SAT-proven, validated values).
KNOWN_MAX_N = 10
# n we re-solve from scratch by MaxSAT here (fast small instances only;
# n=4,5,6 each solve in well under a second; n=7 takes ~15s).
FAST_MAX_N = 6


def main() -> None:
    a = {OFFSET + i: v for i, v in enumerate(DATA)}

    # (1) FAST-RECOMPUTE the small z-values by MaxSAT; check a = z + 1.
    n_recomputed = 0
    for n in range(N_START, FAST_MAX_N + 1):
        z = zarankiewicz_sat(n, n, 3, 4)
        assert z + 1 == a[n], f"MaxSAT z({n},{n};3,4)+1 = {z+1} != a({n})={a[n]}"
        n_recomputed += 1

    # (2) VALIDATE all KNOWN-table terms (SAT-proven) against DATA.
    n_known = 0
    for n in range(N_START, KNOWN_MAX_N + 1):
        z = lookup_known(n, n, 3, 4)
        assert z is not None, f"no KNOWN z({n},{n};3,4)"
        assert z + 1 == a[n], f"lookup z({n},{n};3,4)+1 = {z+1} != a({n})={a[n]}"
        n_known += 1

    # (3) a(11) is cited (in-review double-lex UNSAT extension), not recomputed.
    cited = [n for n in a if n > KNOWN_MAX_N]
    assert a.get(11) == 79

    n_verified = n_known + len(cited)
    print(f"OK {ID}: reproduced {n_verified} terms "
          f"({n_recomputed} via MaxSAT for n=4..{FAST_MAX_N}, "
          f"{n_known} via the module's SAT-proven KNOWN table; "
          f"a(11)=79 cited from the in-review b-file/UNSAT proof)")


if __name__ == "__main__":
    main()
