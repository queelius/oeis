#!/usr/bin/env python3
"""Reproduce and verify A394445: distinct-variable 2-color Rado numbers for
x + y = n*z.

a(n) is the minimum N such that every 2-coloring of {1, ..., N} contains a
monochromatic solution to x + y = n*z with x, y, z pairwise distinct.

Generator: src/rado_numbers.py
  rado_number([1, 1, -n], k_colors=2, distinct=True) -- SAT (CaDiCaL) binary
  search for the least N at which the 2-coloring CNF becomes UNSAT.

Theorem (Towell, 2026): for n >= 8, a(n) = n*(n+3)/2 if n odd,
(n^2 + 2*n + 2)/2 if n even.  (Proved; 500 b-file terms.)

Strategy (kept under ~60s):
  * FAST-RECOMPUTE a(n) by SAT for n = 1..20 (small universes; ~2s total),
    incl. the non-monotone a(7)=49 > a(8)=41.  Check against the b-file.
  * VERIFY every b-file term with n >= 8 against the proven closed form.
  * CONFIRM the published 43-term DATA equals the b-file prefix.

Run as:  python3 reproduce.py   (reads b394445.txt from this directory)
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.rado_numbers import rado_number  # noqa: E402

ID = "A394445"
OFFSET = 1
HERE = os.path.dirname(__file__)
BFILE = os.path.join(HERE, "b394445.txt")

# Published DATA (43 terms, a(1..43)) from the OEIS submission.
DATA = [9, 9, 15, 20, 25, 31, 49, 41, 54, 61, 77, 85, 104, 113, 135, 145, 170,
        181, 209, 221, 252, 265, 299, 313, 350, 365, 405, 421, 464, 481, 527,
        545, 594, 613, 665, 685, 740, 761, 819, 841, 902, 925, 989]

FAST_MAX_N = 20  # recompute a(1..20) from scratch by SAT (fast)


def closed_form(n: int) -> int:
    """Proven closed form for n >= 8."""
    if n % 2 == 1:
        return n * (n + 3) // 2
    return (n * n + 2 * n + 2) // 2


def read_bfile(path: str) -> dict[int, int]:
    """Parse an OEIS b-file ('<n> <a(n)>' per line) into {n: a(n)}."""
    out: dict[int, int] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            out[int(parts[0])] = int(parts[1])
    return out


def main() -> None:
    bf = read_bfile(BFILE)

    # (0) Published DATA must equal the b-file prefix a(1..43).
    for i, v in enumerate(DATA):
        n = OFFSET + i
        assert bf[n] == v, f"b-file a({n})={bf[n]} != published DATA {v}"

    # (1) FAST-RECOMPUTE a(1..20) by SAT; check against the b-file.
    n_recomputed = 0
    for n in range(OFFSET, FAST_MAX_N + 1):
        val, _ = rado_number([1, 1, -n], k_colors=2, max_n=400, distinct=True)
        assert val == bf[n], f"SAT a({n})={val} != b-file {bf[n]}"
        n_recomputed += 1

    # (2) VERIFY every b-file term with n >= 8 against the proven closed form.
    n_formula = 0
    for n, v in sorted(bf.items()):
        if n >= 8:
            assert closed_form(n) == v, \
                f"closed form a({n})={closed_form(n)} != b-file {v}"
            n_formula += 1

    print(f"OK {ID}: reproduced {len(bf)} terms "
          f"({n_recomputed} recomputed by SAT for n=1..{FAST_MAX_N}; "
          f"{n_formula} b-file terms (n>=8) verified against the proven "
          f"closed form; published 43-term DATA matches the b-file)")


if __name__ == "__main__":
    main()
