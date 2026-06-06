#!/usr/bin/env python3
"""Reproduce S(Z/nZ, 3): maximum 3-colorable sum-free subset of Z/nZ.

a(n) = the largest size of a subset A of Z/nZ that can be partitioned into
three sum-free sets (equivalently 3-colored so that no color class contains
a, b, c with a + b == c mod n).

Run standalone from this directory:

    python3 reproduce.py

Recomputes a(n) for n = 2..15 from the value-generating code copied into
``src/schur_backtrack.py`` (the function ``backtrack_schur(orders, k=3)``,
a backtracking + forward-checking solver on the fixed Schur-triple constraint
graph), asserts the result matches the known data, and prints a success line.

The brute-force ``schur_number(orders, k=3)`` in ``src/schur_groups.py``
returns identical values but enumerates all 3^|subset| colorings and becomes
intractable past n ~ 11; ``backtrack_schur`` reaches n = 15 in well under a
second. We cross-check the two solvers on the tractable range n = 2..10.

Runtime: < 1 s.
"""

import sys
import time
from pathlib import Path

# Make the oeis repo root importable so ``from src.<module>`` resolves.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.schur_backtrack import backtrack_schur  # noqa: E402
from src.schur_groups import cyclic_group, schur_number  # noqa: E402

# Known data, offset 2 (n = 2, 3, ..., 15).
OFFSET = 2
KNOWN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13]


def compute(n_max: int = 15):
    """Return [a(2), a(3), ..., a(n_max)] via the fast backtracking solver."""
    return [backtrack_schur(cyclic_group(n), k=3) for n in range(OFFSET, n_max + 1)]


def main() -> int:
    t0 = time.time()
    n_max = OFFSET + len(KNOWN) - 1  # = 15
    terms = compute(n_max)
    assert terms == KNOWN, (
        f"mismatch:\n  computed = {terms}\n  known    = {KNOWN}"
    )

    # Independent cross-check: the brute-force k=3 solver in schur_groups
    # must agree on the range where it is still tractable (n <= 10).
    brute = [schur_number(cyclic_group(n), k=3) for n in range(OFFSET, 11)]
    assert brute == KNOWN[:len(brute)], (
        f"brute/backtrack disagree:\n  brute     = {brute}\n  "
        f"backtrack = {KNOWN[:len(brute)]}"
    )

    elapsed = time.time() - t0
    print(f"OK S-Zn-3: reproduced {len(terms)} terms "
          f"(n={OFFSET}..{n_max}) in {elapsed:.1f}s "
          f"(cross-checked brute force n=2..10)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
