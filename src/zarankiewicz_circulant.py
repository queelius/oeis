"""Diagonal-circulant lower bound for the Zarankiewicz number z(m, n; s, t).

A "diagonal-circulant" m x n binary matrix is one of the form
    M[i, j] = f((i + j) mod n)
for some indicator f: Z_n -> {0, 1}. Each row has the same weight w(f);
each column j has weight w(f) - f[(j - 1) mod n] (with j > 0) or
w(f) (for j = 0). Total ones approximately m * w(f).

A 3x4 submatrix (rows R, columns C) is all-ones iff
    f[(r + c) mod n] = 1 for all (r, c) in R x C,
i.e., R + C (Minkowski sum mod n) is a subset of supp(f).

So we want max |S| (S = supp(f) subset of Z_n) such that no Minkowski
sum R + C with |R| = s, |C| = t lies entirely in S. Equivalently:
    S^c is a hitting set for the family {R + C : R in C_s, C in C_t}.

For n = 11 and z(10, 11; 3, 4) the search is over 2^11 = 2048 subsets,
fast.

The resulting witness gives a constructive lower bound for z(m, n; s, t).
For z(10, 11; 3, 4) the published lower bound is at least 66 (from the
catalog); we check whether circulants reach or exceed that.
"""

from __future__ import annotations

import argparse
import json
import time
from itertools import combinations
from pathlib import Path


def forbidden_minkowski_sums(m: int, n: int, s: int, t: int) -> list[frozenset[int]]:
    """All distinct Minkowski sums R + C (mod n) with R in C(m, s), C in C(n, t)."""
    sums: set[frozenset[int]] = set()
    for R in combinations(range(m), s):
        for C in combinations(range(n), t):
            D = frozenset((r + c) % n for r in R for c in C)
            sums.add(D)
    return sorted(sums, key=lambda d: (len(d), tuple(sorted(d))))


def best_diag_circulant(m: int, n: int, s: int, t: int) -> dict:
    """Find max-weight f: Z_n -> {0, 1} avoiding all-ones s x t submatrix."""
    t0 = time.time()
    forbidden = forbidden_minkowski_sums(m, n, s, t)
    # For each f, supp(f) must not contain any D.
    # Enumerate 2^n subsets in decreasing weight order; first valid is the max.
    best_w = -1
    best_S: list[int] = []
    for w_target in range(n, -1, -1):
        for combo in combinations(range(n), w_target):
            S = frozenset(combo)
            if all(not D.issubset(S) for D in forbidden):
                best_w = w_target
                best_S = sorted(S)
                elapsed = time.time() - t0
                return {
                    "m": m, "n": n, "s": s, "t": t,
                    "max_w_f": best_w,
                    "S": best_S,
                    "row_weight": best_w,
                    "total_ones_lower_bound": m * best_w,
                    "forbidden_count": len(forbidden),
                    "elapsed_s": elapsed,
                }
        # If no S of size w_target works, try smaller.
    # Should always have at least w=0 (S empty avoids all D).
    elapsed = time.time() - t0
    return {
        "m": m, "n": n, "s": s, "t": t,
        "max_w_f": 0,
        "S": [],
        "row_weight": 0,
        "total_ones_lower_bound": 0,
        "forbidden_count": len(forbidden),
        "elapsed_s": elapsed,
    }


def materialize_matrix(m: int, n: int, S: list[int]) -> list[list[int]]:
    """Build the m x n diagonal-circulant matrix from supp(f) = S."""
    Sset = set(S)
    return [[1 if ((i + j) % n) in Sset else 0 for j in range(n)] for i in range(m)]


def verify_no_st_submatrix(M: list[list[int]], s: int, t: int) -> bool:
    """True iff no s x t all-ones submatrix exists in M."""
    m, n = len(M), len(M[0])
    for R in combinations(range(m), s):
        for C in combinations(range(n), t):
            if all(M[r][c] == 1 for r in R for c in C):
                return False
    return True


def count_ones(M: list[list[int]]) -> int:
    return sum(sum(row) for row in M)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("m", type=int, help="rows")
    ap.add_argument("n", type=int, help="columns (modulus of cyclic action)")
    ap.add_argument("s", type=int, help="forbidden submatrix rows")
    ap.add_argument("t", type=int, help="forbidden submatrix columns")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    print(f"Diagonal-circulant z({args.m}, {args.n}; {args.s}, {args.t}) lower bound",
          flush=True)
    result = best_diag_circulant(args.m, args.n, args.s, args.t)
    print(f"  forbidden Minkowski sums: {result['forbidden_count']}", flush=True)
    print(f"  max row-weight f: {result['max_w_f']}", flush=True)
    print(f"  S = supp(f): {result['S']}", flush=True)
    print(f"  total ones: {result['total_ones_lower_bound']} "
          f"(= m * w = {args.m} * {result['max_w_f']})",
          flush=True)
    print(f"  elapsed: {result['elapsed_s']:.3f}s", flush=True)

    if result["max_w_f"] > 0:
        M = materialize_matrix(args.m, args.n, result["S"])
        ok = verify_no_st_submatrix(M, args.s, args.t)
        ones = count_ones(M)
        print(f"  direct verification: no {args.s}x{args.t} all-ones submatrix: {ok}",
              flush=True)
        print(f"  direct ones count: {ones}", flush=True)
        result["matrix"] = M
        result["verified"] = ok
        result["ones_count"] = ones

    if args.out:
        out_path = Path(args.out)
    else:
        out_dir = Path("discoveries/extremal-graph-theory/zarankiewicz/circulant")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"z_{args.m}_{args.n}_{args.s}_{args.t}_circulant.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(f"  wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
