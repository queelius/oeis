#!/usr/bin/env python3
"""Reproduce S(Z/nZ, 2): maximum 2-colorable sum-free subset of Z/nZ.

a(n) = the largest size of a subset A of Z/nZ that can be partitioned into
two sum-free sets (equivalently 2-colored so that neither color class
contains a, b, c with a + b == c mod n).

Run standalone from this directory:

    python3 reproduce.py

Recomputes a(n) for n = 2..20 from the value-generating code copied into
``src/schur_groups.py`` (the function ``schur_number(orders, k=2)``), asserts
the result matches the known data, and prints a success line.

Runtime: ~20 s (the n = 19 exhaustive search dominates).
"""

import sys
import time
from pathlib import Path

# Make the oeis repo root importable so ``from src.schur_groups`` resolves.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.schur_groups import cyclic_group, schur_number  # noqa: E402

# Known data, offset 2 (n = 2, 3, ..., 20).
OFFSET = 2
KNOWN = [1, 2, 3, 4, 4, 4, 6, 6, 8, 8, 9, 8, 9, 12, 12, 12, 12, 12, 16]


def compute(n_max: int = 20):
    """Return [a(2), a(3), ..., a(n_max)]."""
    return [schur_number(cyclic_group(n), k=2) for n in range(OFFSET, n_max + 1)]


def main() -> int:
    t0 = time.time()
    n_max = OFFSET + len(KNOWN) - 1  # = 20
    terms = compute(n_max)
    assert terms == KNOWN, (
        f"mismatch:\n  computed = {terms}\n  known    = {KNOWN}"
    )
    elapsed = time.time() - t0
    print(f"OK S-Zn-2: reproduced {len(terms)} terms "
          f"(n={OFFSET}..{n_max}) in {elapsed:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
