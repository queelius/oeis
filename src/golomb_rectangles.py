"""Exact computation of Golomb-rectangle numbers G(m, n).

A *Golomb rectangle* is a placement of dots in an m x n grid such that all
pairwise difference vectors between dots are distinct -- a two-dimensional
Sidon set (B_2 set) in the grid.  Equivalently (Costas-array survey,
Drakakis & Iorgov / Golomb & Taylor, arXiv:1102.5727, Definition 12 and
Problem 23, after J. P. Robinson, "Golomb rectangles", IEEE Trans. IT, 1985):
overlay two copies of the 0/1 dot pattern and slide one by any nonzero
offset; the dots overlap in at most one position.  The autocorrelation range
is then {N, 1, 0} where N is the number of dots.

G(m, n) is the MAXIMUM number of dots placeable in an m x n grid with the
distinct-difference property.  This is the extremal ("maximal N") side of
Problem 23 of the Costas survey, which states that "virtually nothing is
known about the mathematical properties of these objects."  There is no OEIS
array for G(m, n); only the one-dimensional optimal Golomb ruler A003022.

The difference-vector condition, made precise
--------------------------------------------
Let the dots be points p_1, ..., p_k in Z^2.  Write the ordered difference of
an ordered pair (p_i, p_j), i != j, as p_i - p_j.  The set is a Golomb
rectangle iff all these ordered differences are distinct.  Since
p_i - p_j = -(p_j - p_i) and a nonzero vector d != -d, this is the same as
requiring all C(k, 2) UNORDERED-pair differences to be distinct up to sign;
it is also the same as the autocorrelation condition above.  A repeated
difference is witnessed by four dots A, B, C, D with
    B - A = D - C = d != 0   and   (A, B) != (C, D),
so forbidding the simultaneous presence of every such "bad quadruple"
exactly enforces the property.

Methods
-------
  - brute_force_max(m, n): independent ground truth by enumerating point sets
    with a growth/pruning search (no SAT).  Used to validate the SAT path.
  - sat_can_place(m, n, k): decision oracle "can >= k dots be placed?" via a
    CNF that forbids every bad quadruple, plus a cardinality bound, plus
    optional grid-dihedral lex-leader symmetry breaking.
  - golomb_rectangle(m, n): G(m, n) by SAT decision + climb, returning the
    value, an optimal witness, and a proof flag (SAT witness at G and UNSAT
    at G + 1, hence proven optimal).

Symmetry breaking
-----------------
The maximum is invariant under the dihedral symmetry of the grid: for m != n
the Klein four-group {identity, horizontal flip, vertical flip, 180 rotation}
(order 4); for the square m == n the full dihedral group D4 (order 8), which
additionally permits transpose and the two diagonal reflections.

IMPORTANT (why double-lex is WRONG here): the distinct-difference property is
NOT invariant under arbitrary row/column permutations -- only under these
RIGID grid motions.  (E.g. a 4x4 Golomb rectangle with 6 dots loses the
property after swapping just rows 0 and 1.)  This is the crucial difference
from the Zarankiewicz "no all-ones submatrix" property, which IS invariant
under every row and column permutation and therefore admits the double-lex
break of src/zarankiewicz_symbreak.py.  Applying double-lex to the Sidon
property is UNSOUND: it can discard an entire orbit whose lexicographic leader
is a non-Golomb permutation, undercounting G (observed: it wrongly returns
G(4,4)=5 instead of 6).

We instead break exactly the grid's dihedral symmetry: we require the chosen
0/1 grid M to be lexicographically smallest (in row-major flattening) among
its dihedral images sigma(M).  For each non-identity symmetry sigma in the
group (D2 or D4), we add the lex constraint  M <=_lex sigma(M).  Every orbit
under the dihedral group has a unique lex-least member, so this is sound:
a witness exists under the constraint iff one exists at all.

Soundness summary:
  * SAT at k    <=>  G(m, n) >= k
  * UNSAT at k  <=>  G(m, n) <  k
both with and without symmetry breaking.
"""

from __future__ import annotations

import argparse
from itertools import combinations

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver


# ---------------------------------------------------------------------------
# Known values (for validation only -- NOT used by the solver).
#
# G(m, n) is symmetric: G(m, n) = G(n, m).  Trivial corners:
#   G(1, n): a single row is a 1D Golomb ruler living in length-n; the maximum
#            number of marks is A003022's inverse, i.e. the largest k whose
#            optimal ruler length <= n - 1.  For n = 1..6 this is
#            1, 2, 3, 3, 4, 4 (rulers {0},{0,1},{0,1,3},{0,1,3} in len 3,
#            {0,1,4,6} needs length 6 so n=7... see below).
# We seed only the handful of values we cross-check by brute force in tests;
# the solver derives everything from scratch.
# ---------------------------------------------------------------------------

# Small verified values G(m, n) for 2 <= m <= n, cross-checked by brute force
# (asserted in reproduce.py and tests; kept here for reference).  These were
# established by independent brute-force enumeration, NOT assumed.
KNOWN_G = {
    (2, 2): 3,
    (2, 3): 4,
    (2, 4): 4,
    (3, 3): 5,
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def cells(m: int, n: int):
    """All grid cells (r, c), 0 <= r < m, 0 <= c < n, row-major."""
    return [(r, c) for r in range(m) for c in range(n)]


def _has_distinct_differences(points) -> bool:
    """True iff all ordered pairwise difference vectors among `points` are
    distinct (equivalently all unordered differences distinct up to sign)."""
    seen = set()
    pts = list(points)
    for i in range(len(pts)):
        ri, ci = pts[i]
        for j in range(len(pts)):
            if i == j:
                continue
            d = (ri - pts[j][0], ci - pts[j][1])
            if d in seen:
                return False
            seen.add(d)
    return True


# ---------------------------------------------------------------------------
# Brute force (independent ground truth)
# ---------------------------------------------------------------------------

def brute_force_max(m: int, n: int):
    """Return (G, witness) by an exhaustive growth search over dot sets.

    Builds Golomb-rectangle dot sets incrementally in row-major cell order,
    pruning the moment the distinct-difference property breaks.  This is a
    pure combinatorial enumeration with NO SAT, used to validate the SAT path.
    Returns the maximum size found and one realizing point set.
    """
    grid = cells(m, n)
    best = {"size": 0, "pts": []}

    def diffs_of(pts):
        s = set()
        for a in range(len(pts)):
            for b in range(len(pts)):
                if a == b:
                    continue
                s.add((pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]))
        return s

    def extend(start_idx, chosen, used_diffs):
        if len(chosen) > best["size"]:
            best["size"] = len(chosen)
            best["pts"] = list(chosen)
        # Bound: even taking every remaining cell cannot beat best.
        if len(chosen) + (len(grid) - start_idx) <= best["size"]:
            return
        for idx in range(start_idx, len(grid)):
            p = grid[idx]
            new = set()
            ok = True
            for q in chosen:
                d1 = (p[0] - q[0], p[1] - q[1])
                d2 = (q[0] - p[0], q[1] - p[1])
                if d1 in used_diffs or d2 in used_diffs or d1 in new or d2 in new:
                    ok = False
                    break
                new.add(d1)
                new.add(d2)
            if ok:
                chosen.append(p)
                used_diffs |= new
                extend(idx + 1, chosen, used_diffs)
                used_diffs -= new
                chosen.pop()

    extend(0, [], set())
    return best["size"], best["pts"]


def brute_force_count_at(m: int, n: int, k: int) -> int:
    """Number of k-dot Golomb rectangles in the m x n grid (no symmetry
    reduction; counts labelled placements).  Used only for tiny sanity checks
    and C(m, n, N) spot-values."""
    grid = cells(m, n)
    count = 0
    for combo in combinations(grid, k):
        if _has_distinct_differences(combo):
            count += 1
    return count


# ---------------------------------------------------------------------------
# SAT encoding
# ---------------------------------------------------------------------------

def _bad_quadruples(m: int, n: int):
    """Yield frozensets of cells {A, B, C, D} whose simultaneous presence
    forces a repeated difference: B - A = D - C = d != 0, (A,B) != (C,D).

    We enumerate by difference vector d: for each nonzero d realizable in the
    grid, list all ordered pairs (P, P+d) inside the grid, then take every
    unordered pair of two DISTINCT such pairs.  The involved cells form a set
    of size 4 (generic) or size 3 (when P1+d = P2, i.e. the colinear chain
    P1, P1+d, P1+2d realizes d twice through a shared middle dot); size 2 is
    impossible since it would require 2d = 0.  That set must not be all-dots.
    """
    grid_set = set(cells(m, n))
    # Distinct nonzero difference vectors with dr >= 0 (and dc > 0 if dr == 0);
    # each unordered direction once -- the pair list below already covers both
    # orientations because we iterate ordered (P, P+d).
    diffs = set()
    for r in range(-(m - 1), m):
        for c in range(-(n - 1), n):
            if (r, c) == (0, 0):
                continue
            if (r > 0) or (r == 0 and c > 0):
                diffs.add((r, c))

    seen = set()
    for (dr, dc) in diffs:
        pairs = []
        for (r, c) in grid_set:
            q = (r + dr, c + dc)
            if q in grid_set:
                pairs.append(((r, c), q))
        for (p1, p2) in combinations(pairs, 2):
            quad = frozenset([p1[0], p1[1], p2[0], p2[1]])
            if quad not in seen:
                seen.add(quad)
                yield quad


def build_encoding(m: int, n: int, k: int, symmetry: bool = True):
    """Build the CNF deciding "can >= k dots be placed?".

    Returns (clauses, var_of, pool) where var_of[(r, c)] is the cell variable.
    Variables: one per cell, then auxiliary cardinality/lex vars via IDPool.
    """
    pool = IDPool()
    var_of = {}
    for (r, c) in cells(m, n):
        var_of[(r, c)] = pool.id(("x", r, c))

    clauses: list[list[int]] = []

    # (1) Distinct differences: forbid every bad quadruple.
    for quad in _bad_quadruples(m, n):
        clauses.append([-var_of[cell] for cell in quad])

    # (2) Cardinality: at least k dots.
    cell_vars = [var_of[(r, c)] for (r, c) in cells(m, n)]
    card = CardEnc.atleast(lits=cell_vars, bound=k, vpool=pool,
                           encoding=EncType.seqcounter)
    clauses.extend(card.clauses)

    # (3) Symmetry breaking: M is lex-least among its dihedral images.
    if symmetry:
        clauses.extend(_dihedral_symmetry_clauses(m, n, var_of, pool))

    return clauses, var_of, pool


def dihedral_images(m: int, n: int):
    """Yield (name, mapping) for each symmetry of the m x n grid that
    preserves the distinct-difference property.

    `mapping` sends a source cell (r, c) -> image cell.  For m != n the group
    is {id, flip-rows, flip-cols, rot180} (order 4); for the square m == n it
    also includes {transpose, anti-transpose, rot90, rot270} (order 8).  All
    images stay within the SAME m x n grid (transpose maps n x n -> n x n).
    """
    base = [
        ("id", lambda r, c: (r, c)),
        ("flipR", lambda r, c: (m - 1 - r, c)),
        ("flipC", lambda r, c: (r, n - 1 - c)),
        ("rot180", lambda r, c: (m - 1 - r, n - 1 - c)),
    ]
    if m == n:
        base += [
            ("transpose", lambda r, c: (c, r)),
            ("antitrans", lambda r, c: (n - 1 - c, n - 1 - r)),
            ("rot90", lambda r, c: (c, n - 1 - r)),
            ("rot270", lambda r, c: (n - 1 - c, r)),
        ]
    return base


def _dihedral_symmetry_clauses(m: int, n: int, var_of, pool):
    """Require M (row-major flattening) to be lex-least over its dihedral
    orbit: for every non-identity symmetry sigma, add  M <=_lex sigma(M).

    Sound for the Sidon property because the dihedral group is exactly the set
    of grid motions preserving distinct differences, and every orbit has a
    unique lex-least representative.
    """
    order = cells(m, n)  # fixed row-major cell order for flattening
    a = [var_of[cell] for cell in order]
    clauses = []
    for name, mapping in dihedral_images(m, n):
        if name == "id":
            continue
        b = [var_of[mapping(r, c)] for (r, c) in order]
        clauses.extend(_lex_leq_clauses(a, b, pool, tag=("sym", name)))
    return clauses


def _lex_leq_clauses(a, b, pool, tag):
    """CNF for a <=_lex b between equal-length variable lists (prefix-equal
    auxiliary encoding; cf. src/zarankiewicz_symbreak.lex_leq_clauses).

    `tag` is a hashable identifier UNIQUE to this lex constraint; it keys the
    fresh auxiliary variables in the shared IDPool.  (Using id() of the python
    lists is unsafe -- transient lists can share an id, which silently aliases
    auxiliary vars across constraints and corrupts the encoding.)
    """
    assert len(a) == len(b)
    clauses = []
    prefix = None  # None == constant-true (empty equal prefix)
    kk = len(a)
    for j in range(kk):
        aj, bj = a[j], b[j]
        if prefix is None:
            clauses.append([-aj, bj])
        else:
            clauses.append([-prefix, -aj, bj])
        if j == kk - 1:
            break
        eqj = pool.id(("eq", tag, j))
        clauses.append([-eqj, -aj, bj])
        clauses.append([-eqj, -bj, aj])
        clauses.append([aj, bj, eqj])
        clauses.append([-aj, -bj, eqj])
        if prefix is None:
            prefix = eqj
        else:
            pj1 = pool.id(("pfx", tag, j))
            clauses.append([-pj1, prefix])
            clauses.append([-pj1, eqj])
            clauses.append([-prefix, -eqj, pj1])
            prefix = pj1
    return clauses


def sat_can_place(m: int, n: int, k: int, symmetry: bool = True,
                  solver_name: str = "g3"):
    """Decide whether >= k dots fit.  Returns (sat: bool, witness or None)."""
    if k <= 1:
        # 0 or 1 dot is always a valid Golomb rectangle (no pairs).
        pts = []
        if k == 1:
            pts = [(0, 0)]
        return True, pts
    clauses, var_of, _pool = build_encoding(m, n, k, symmetry=symmetry)
    with Solver(name=solver_name, bootstrap_with=clauses) as solver:
        sat = solver.solve()
        if not sat:
            return False, None
        model = set(v for v in solver.get_model() if v > 0)
        witness = [(r, c) for (r, c) in cells(m, n) if var_of[(r, c)] in model]
        return True, witness


# ---------------------------------------------------------------------------
# G(m, n) by SAT decision + climb
# ---------------------------------------------------------------------------

def golomb_rectangle(m: int, n: int, symmetry: bool = True,
                     solver_name: str = "g3", lower: int = 1):
    """Compute G(m, n) exactly via SAT decision + climb.

    Climbs k = max(lower, 1), lower+1, ... while SAT, keeping the last witness.
    G is the largest k with SAT; optimality is *proven* because k = G + 1 is
    UNSAT.  Returns a dict: {m, n, G, witness, proven, symmetry}.

    `lower` lets a known lower bound (e.g. from brute force or a construction)
    skip the easy SAT calls; it does not affect the proven flag.
    """
    if m > n:
        m, n = n, m  # G is symmetric; canonicalize m <= n
    best_k = 0
    best_witness: list = []
    proven = False
    k = max(1, lower)
    # Make sure we have a valid baseline at k-1 (so best reflects a witness).
    if k > 1:
        sat0, w0 = sat_can_place(m, n, k - 1, symmetry=symmetry,
                                 solver_name=solver_name)
        if sat0:
            best_k = k - 1
            best_witness = w0
    while True:
        sat, witness = sat_can_place(m, n, k, symmetry=symmetry,
                                     solver_name=solver_name)
        if sat:
            best_k = k
            best_witness = witness
            k += 1
        else:
            proven = True  # k is UNSAT, so G = k - 1 is optimal
            break
        # Hard cap: cannot exceed total cells.
        if k > m * n:
            proven = True
            break
    return {
        "m": m, "n": n, "G": best_k, "witness": sorted(best_witness),
        "proven": proven, "symmetry": symmetry,
    }


def verify_witness(m: int, n: int, points) -> bool:
    """Independently verify a claimed Golomb-rectangle witness: all points in
    grid, distinct, and all difference vectors distinct."""
    pts = list(points)
    if any(not (0 <= r < m and 0 <= c < n) for (r, c) in pts):
        return False
    if len(set(map(tuple, pts))) != len(pts):
        return False
    return _has_distinct_differences(pts)


# ---------------------------------------------------------------------------
# Table builder
# ---------------------------------------------------------------------------

def build_table(max_m: int, max_n: int, symmetry: bool = True,
                solver_name: str = "g3", verbose: bool = False):
    """Compute G(m, n) for 1 <= m <= n <= max, with m <= max_m, n <= max_n.

    Returns a dict {(m, n): result_dict}.  Uses the brute-force lower bound to
    warm-start the climb only for the smallest cases; otherwise climbs from 1.
    """
    table = {}
    for m in range(1, max_m + 1):
        for n in range(m, max_n + 1):
            res = golomb_rectangle(m, n, symmetry=symmetry,
                                   solver_name=solver_name)
            assert verify_witness(m, n, res["witness"]), \
                f"witness failed verification at ({m},{n})"
            table[(m, n)] = res
            if verbose:
                flag = "proven" if res["proven"] else "lower-bound"
                print(f"G({m},{n}) = {res['G']:2d}  [{flag}]  "
                      f"witness={res['witness']}")
    return table


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli():
    ap = argparse.ArgumentParser(description="Golomb rectangle G(m,n) solver")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("one", help="compute a single G(m,n)")
    p1.add_argument("m", type=int)
    p1.add_argument("n", type=int)
    p1.add_argument("--no-sym", action="store_true")
    p1.add_argument("--solver", default="g3")

    p2 = sub.add_parser("table", help="compute a table up to max_m x max_n")
    p2.add_argument("max_m", type=int)
    p2.add_argument("max_n", type=int)
    p2.add_argument("--no-sym", action="store_true")
    p2.add_argument("--solver", default="g3")

    p3 = sub.add_parser("brute", help="brute-force G(m,n) (ground truth)")
    p3.add_argument("m", type=int)
    p3.add_argument("n", type=int)

    args = ap.parse_args()
    if args.cmd == "one":
        res = golomb_rectangle(args.m, args.n, symmetry=not args.no_sym,
                               solver_name=args.solver)
        ok = verify_witness(res["m"], res["n"], res["witness"])
        flag = "PROVEN optimal" if res["proven"] else "lower bound only"
        print(f"G({res['m']},{res['n']}) = {res['G']}  ({flag}); "
              f"witness verified: {ok}")
        print(f"  witness: {res['witness']}")
    elif args.cmd == "table":
        table = build_table(args.max_m, args.max_n,
                            symmetry=not args.no_sym,
                            solver_name=args.solver, verbose=True)
        print(f"\n{len(table)} entries computed.")
    elif args.cmd == "brute":
        g, w = brute_force_max(args.m, args.n)
        print(f"G({args.m},{args.n}) = {g} (brute force); witness={sorted(w)}")


if __name__ == "__main__":
    _cli()
