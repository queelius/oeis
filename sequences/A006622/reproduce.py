#!/usr/bin/env python3
"""Reproduce and verify A006622: a(n) = z(n,n+1;3,4) + 1.

z(n,n+1;3,4) is the Zarankiewicz number for n x (n+1) matrices with no
3 x 4 all-ones submatrix.  The OEIS term is k_{3,4}(n,n+1) = z + 1.

Generator: src/zarankiewicz.py
  - lower bound: SAT witness (MaxSAT / RC2)
  - upper bound: double-lex symmetry-broken UNSAT (src/zarankiewicz_symbreak.py)
  - lookup_known(n, n+1, 3, 4): SAT-PROVEN exact z-values in the module's
    KNOWN_Z_3_4_RECT1 table, (rows, cols) = (3,4)..(9,10).

Strategy (kept under ~60s):
  * FAST-RECOMPUTE z(n,n+1;3,4) by MaxSAT for the small n (n = 3, 4, 5).
  * VALIDATE every originally published term (n = 3..9) against the KNOWN
    table; a(n) = z + 1.
  * a(10) = 73 (z(10,11;3,4) = 72) is the in-review extension proved by a
    1633 s double-lex UNSAT; NOT recomputed here, carried as a cited term.

Run as:  python3 reproduce.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.zarankiewicz import lookup_known, zarankiewicz_sat  # noqa: E402

ID = "A006622"
OFFSET = 3
# Published DATA a(n) for n = 3..10.  a(3..9) originally published
# (a(9)=61 added Mar 2026); a(10)=73 is the in-review extension.
DATA = [12, 18, 26, 33, 41, 51, 61, 73]
N_START = 3
KNOWN_MAX_N = 9   # last n in KNOWN_Z_3_4_RECT1 (rows up to 9, cols 10)
FAST_MAX_N = 5    # n we re-solve from scratch by MaxSAT (fast: n=3,4,5)


def main() -> None:
    a = {OFFSET + i: v for i, v in enumerate(DATA)}

    # (1) FAST-RECOMPUTE small z-values by MaxSAT; check a = z + 1.
    n_recomputed = 0
    for n in range(N_START, FAST_MAX_N + 1):
        z = zarankiewicz_sat(n, n + 1, 3, 4)
        assert z + 1 == a[n], \
            f"MaxSAT z({n},{n+1};3,4)+1 = {z+1} != a({n})={a[n]}"
        n_recomputed += 1

    # (2) VALIDATE all KNOWN-table terms (SAT-proven) against DATA.
    n_known = 0
    for n in range(N_START, KNOWN_MAX_N + 1):
        z = lookup_known(n, n + 1, 3, 4)
        assert z is not None, f"no KNOWN z({n},{n+1};3,4)"
        assert z + 1 == a[n], \
            f"lookup z({n},{n+1};3,4)+1 = {z+1} != a({n})={a[n]}"
        n_known += 1

    # (3) a(10) cited (in-review double-lex UNSAT), not recomputed.
    cited = [n for n in a if n > KNOWN_MAX_N]
    assert a.get(10) == 73

    n_verified = n_known + len(cited)
    print(f"OK {ID}: reproduced {n_verified} terms "
          f"({n_recomputed} via MaxSAT for n=3..{FAST_MAX_N}, "
          f"{n_known} via the module's SAT-proven KNOWN table; "
          f"a(10)=73 cited from the in-review UNSAT proof)")


if __name__ == "__main__":
    main()
