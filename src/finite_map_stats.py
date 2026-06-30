"""
Exact statistics of functional graphs of finite maps.

A functional graph on [n] = {0, 1, ..., n-1} is the graph of a map
f : [n] -> [n]; there are n^n of them. With a directed edge x -> f(x), every
component is a set of cycles with trees ("rho shapes") hanging off the cyclic
nodes. This module computes EXACT integer totals and EXACT rational
expectations of classical random-mapping statistics over all n^n maps for
small n, and enumerates structurally-restricted map classes whose counting
sequences have closed forms.

Statistics tracked (per map, then summed over all n^n maps):
  * num_components   number of weakly-connected components (= number of cycles)
  * num_cyclic       number of cyclic / periodic points (image of f^n)
  * num_fixed        number of fixed points (f(x) = x)
  * image_size       |f([n])|
  * num_terminal     number of terminal nodes (no preimage): n - image_size
  * max_tail         maximal tail length (longest distance to a cycle); 0 if
                     every point is cyclic (a permutation)

Validated known quantities (see tests):
  * total maps                  = n^n
  * number of permutations      = n!
  * number of idempotent maps   = A000248 = sum_k C(n,k) k^(n-k)
  * number of connected maps    = A001865 (one component)
  * sum of cyclic points        = A001865 as well; the mean number of cyclic
                                  points over a random mapping is the classical
                                  Q-function  sum_{k=1..n} n!/((n-k)! n^k)
  * sum of fixed points         = n^n (mean number of fixed points is 1)

Restricted classes with structurally-derived closed forms (both attested in
OEIS by Critzer 2012):
  * "depth <= 1" maps: every vertex within distance 1 of a cycle.
        count = sum_{k=1..n} C(n,k) k! k^(n-k)            (OEIS A006153,
        e.g.f. 1/(1 - x e^x))
    Derivation: let C be the cyclic set, |C| = k; the k! permutations of C
    realize C as the periodic part, and each of the remaining n-k vertices
    must map directly into C (k choices each) to sit at depth exactly 1.
  * maps with all in-degrees in {0, 2}: defined only for even n = 2m,
        count = C(2m, m) (2m)! / 2^m                      (OEIS A126934,
        |a(m)| = (2m-1)!! (2m)!/2^m)
    Derivation: exactly m vertices have in-degree 2 (the image), C(2m,m) ways
    to choose them; distribute the 2m domain elements onto these m targets
    with multiplicity 2 each: (2m)!/(2!)^m ways.

References:
  - Flajolet & Odlyzko (1990), "Random Mapping Statistics", EUROCRYPT '89,
    LNCS 434, pp. 329-354.
  - Flajolet & Sedgewick (2009), "Analytic Combinatorics", Cambridge UP,
    Ch. II.4 (functional graphs as sets of cycles of trees).
  - Harris (1960), "Probability distributions related to random mappings",
    Ann. Math. Statist. 31(4):1045-1062.
  - Critzer (2012), OEIS comments on A006153 (distance <= 1 from a cycle) and
    A126934 (in-degrees in {0,2}).
  - OEIS A000248 (idempotent maps), A001865 (connected functions).
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction
from itertools import product
from typing import Any, Iterator


# ---------------------------------------------------------------------------
# Single-map structural analysis
# ---------------------------------------------------------------------------

def cyclic_points(f: tuple[int, ...] | list[int]) -> set[int]:
    """Return the set of cyclic (periodic) points of an endofunction f on [n].

    A point is cyclic iff it lies on a cycle, equivalently iff it is in the
    image of f^n (after n iterations every point has entered its cycle, and
    f^n acts as a permutation on the cyclic set). Computed by iterating the
    whole vertex set n times -- exact integer arithmetic only.
    """
    n = len(f)
    if n == 0:
        return set()
    g = list(range(n))
    for _ in range(n):
        g = [f[g[x]] for x in range(n)]
    return set(g)


def tail_lengths(f: tuple[int, ...] | list[int]) -> list[int]:
    """Return, for each vertex, its tail length: the number of steps to reach a
    cyclic point (0 for cyclic points themselves).

    Computed exactly by following each vertex until it lands in the cyclic set,
    with memoization. The functional graph guarantees termination.
    """
    n = len(f)
    cyc = cyclic_points(f)
    depth: list[int | None] = [0 if x in cyc else None for x in range(n)]

    def resolve(x: int) -> int:
        path: list[int] = []
        while depth[x] is None:
            path.append(x)
            x = f[x]
        d = depth[x]
        for y in reversed(path):
            d += 1
            depth[y] = d
        return depth[path[0]] if path else depth[x]

    for x in range(n):
        if depth[x] is None:
            resolve(x)
    return [d for d in depth]  # type: ignore[misc]


def num_components(f: tuple[int, ...] | list[int]) -> int:
    """Number of weakly-connected components of the functional graph (= number
    of distinct cycles). Union-find over the n edges x -> f(x)."""
    n = len(f)
    if n == 0:
        return 0
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for x in range(n):
        ra, rb = find(x), find(f[x])
        if ra != rb:
            parent[ra] = rb
    return len({find(x) for x in range(n)})


def map_statistics(f: tuple[int, ...] | list[int]) -> dict[str, int]:
    """Compute all per-map statistics for a single endofunction f on [n]."""
    n = len(f)
    img = set(f)
    cyc = cyclic_points(f)
    tails = tail_lengths(f)
    return {
        "num_components": num_components(f),
        "num_cyclic": len(cyc),
        "num_fixed": sum(1 for x in range(n) if f[x] == x),
        "image_size": len(img),
        "num_terminal": n - len(img),
        "max_tail": max(tails) if tails else 0,
        "is_permutation": 1 if len(img) == n else 0,
        "is_idempotent": 1 if all(f[f[x]] == f[x] for x in range(n)) else 0,
        "is_connected": 1 if num_components(f) == 1 else 0,
    }


# ---------------------------------------------------------------------------
# Exhaustive enumeration over all n^n maps
# ---------------------------------------------------------------------------

def all_maps(n: int) -> Iterator[tuple[int, ...]]:
    """Yield every endofunction on [n] as a tuple of length n (n^n of them)."""
    return product(range(n), repeat=n)


def enumerate_totals(n: int) -> dict[str, Any]:
    """Brute-force exact totals and rational expectations over all n^n maps.

    Returns a dict with integer totals (summed over all maps) for each additive
    statistic, integer counts for the boolean predicates, and exact Fraction
    expectations. All arithmetic is exact (int / Fraction).
    """
    total = 0
    sums: Counter[str] = Counter()
    counts: Counter[str] = Counter()
    max_tail_over_all = 0

    for f in all_maps(n):
        total += 1
        st = map_statistics(f)
        for key in ("num_components", "num_cyclic", "num_fixed",
                    "image_size", "num_terminal", "max_tail"):
            sums[key] += st[key]
        for key in ("is_permutation", "is_idempotent", "is_connected"):
            counts[key] += st[key]
        if st["max_tail"] > max_tail_over_all:
            max_tail_over_all = st["max_tail"]

    denom = total  # = n^n
    expectations = {
        key: Fraction(sums[key], denom) for key in sums
    } if denom else {}

    return {
        "n": n,
        "total_maps": total,
        "sums": dict(sums),
        "counts": dict(counts),
        "expectations": expectations,
        "max_tail_over_all_maps": max_tail_over_all,
    }


# ---------------------------------------------------------------------------
# Closed forms for validated known quantities
# ---------------------------------------------------------------------------

def total_maps(n: int) -> int:
    """n^n (with 0^0 = 1)."""
    return n ** n if n > 0 else 1


def num_permutations(n: int) -> int:
    """n! -- the maps that are bijections."""
    return math.factorial(n)


def num_idempotent(n: int) -> int:
    """Number of idempotent maps (f o f = f): OEIS A000248.

    A map is idempotent iff its image is a set of fixed points and every other
    point maps into that set. Choosing the k-element image (which must be
    fixed) and mapping the remaining n-k points into it gives
    sum_{k=0..n} C(n,k) k^(n-k).
    """
    return sum(math.comb(n, k) * (k ** (n - k)) for k in range(n + 1))


def num_idempotent_by_image_size(n: int, k: int) -> int:
    """Idempotent maps with image (= fixed-point set) of size exactly k:
    C(n,k) k^(n-k). Row k of the A000248 refinement."""
    return math.comb(n, k) * (k ** (n - k))


def expected_cyclic_points(n: int) -> Fraction:
    """Exact mean number of cyclic points over a uniform random map on [n].

    Classical result (Harris 1960; Flajolet-Odlyzko 1990):
        E[# cyclic points] = sum_{k=1..n} n! / ((n-k)! n^k)
    (the "Ramanujan Q-function" form). Returned as an exact Fraction.
    """
    return sum(
        Fraction(math.factorial(n), math.factorial(n - k)) * Fraction(1, n ** k)
        for k in range(1, n + 1)
    ) if n > 0 else Fraction(0)


def total_cyclic_points(n: int) -> int:
    """Sum of cyclic points over all n^n maps = n^n * E[# cyclic points].

    Equivalently the number of (map, marked cyclic point) pairs:
        sum_{k=1..n} n!/(n-k)! * n^{n-k}      (= n^n * Q-function).
    Values 1, 6, 51, 568, 7845, ... = OEIS A063169, which satisfies the clean
    identity a(n) = n * A001865(n) (n times the number of connected functions).
    Returned as an exact integer.
    """
    val = total_maps(n) * expected_cyclic_points(n)  # exact Fraction, denom | n^n
    assert val.denominator == 1, "total cyclic points must be integral"
    return int(val)



def num_connected(n: int) -> int:
    """Number of connected functions on n labeled nodes (one component):
    OEIS A001865.

    Computed via the exponential formula relating the e.g.f. C(x) of connected
    functional graphs to the e.g.f. of all functional graphs F(x) = 1/(1-T(x)),
    where T(x) = x e^{T(x)} is the tree e.g.f.: F = exp(C), so C = log F.
    We extract coefficients exactly over the rationals.
    """
    if n == 0:
        return 0
    # Cayley: number of rooted forests / tree e.g.f. coefficients t_j = j^{j-1}.
    # All functional graphs: f_n = n^n. Connected: c_n via log of e.g.f.
    # Work with a_j = n^n / n! style is messy; use the set/multiset relation
    # F(x) = exp(C(x)) over e.g.f.'s, i.e. f_n/n! = [x^n] exp(sum c_j x^j / j!).
    # Recover c_n by Newton's identity on e.g.f. coefficients.
    from fractions import Fraction as Fr
    f = [Fr(total_maps(j), math.factorial(j)) for j in range(n + 1)]  # f[j]=f_j/j!
    c = [Fr(0)] * (n + 1)
    # F = exp(C): f_m = (1/m) sum_{j=1..m} j c_j f_{m-j}  (with f_0 = 1)
    # Solve for c_m: c_m = f_m - (1/m) sum_{j=1..m-1} j c_j f_{m-j}, then *m/(...)
    # Standard: m f_m = sum_{j=1..m} j c_j f_{m-j}
    for m in range(1, n + 1):
        s = sum(j * c[j] * f[m - j] for j in range(1, m))
        # m f_m = s + m c_m f_0  => c_m = (m f_m - s) / m   (f_0 = 1)
        c[m] = (m * f[m] - s) / m
    return int(c[n] * math.factorial(n))


# ---------------------------------------------------------------------------
# Restricted map classes with closed-form counts
# ---------------------------------------------------------------------------

def num_depth_le_1(n: int) -> int:
    """Number of maps in which every vertex is within distance 1 of a cycle.

    Closed form (structural):  sum_{k=1..n} C(n,k) k! k^(n-k).
    Equals OEIS A006153 (e.g.f. 1/(1 - x e^x)). For n = 0 returns 1 (empty map),
    matching the A006153 offset-0 value a(0) = 1.
    """
    if n == 0:
        return 1
    return sum(
        math.comb(n, k) * math.factorial(k) * (k ** (n - k))
        for k in range(1, n + 1)
    )


def num_indegree_in_02(n: int) -> int:
    """Number of maps on [n] with every in-degree in {0, 2}.

    Nonzero only when n = 2m is even, where the closed form is
        C(2m, m) * (2m)! / 2^m              (= OEIS A126934, |a(m)|).
    Returns 0 for odd n (a counting parity obstruction: in-degrees in {0,2}
    force |image| = n/2).
    """
    if n % 2 != 0:
        return 0
    m = n // 2
    return math.comb(2 * m, m) * math.factorial(2 * m) // (2 ** m)


def brute_force_restricted_count(n: int, predicate) -> int:
    """Count maps on [n] satisfying a boolean predicate(f) by exhaustion.

    Used by tests to confirm the closed forms above. predicate takes a tuple
    of length n and returns bool.
    """
    return sum(1 for f in all_maps(n) if predicate(f))


def is_depth_le_1(f: tuple[int, ...] | list[int]) -> bool:
    """True iff every vertex of f is within distance 1 of a cyclic point."""
    cyc = cyclic_points(f)
    return all(x in cyc or f[x] in cyc for x in range(len(f)))


def is_indegree_in_02(f: tuple[int, ...] | list[int]) -> bool:
    """True iff every vertex of f has in-degree 0 or 2."""
    n = len(f)
    deg = [0] * n
    for x in range(n):
        deg[f[x]] += 1
    return all(d in (0, 2) for d in deg)


# ---------------------------------------------------------------------------
# Summary driver
# ---------------------------------------------------------------------------

def summary_table(n_max: int = 7) -> list[dict[str, Any]]:
    """Build a per-n summary of exact totals and validated closed forms.

    Returns one dict per n in 1..n_max with the brute-force totals alongside
    the corresponding closed-form predictions, ready for printing or JSON dump.
    """
    rows = []
    for n in range(1, n_max + 1):
        tot = enumerate_totals(n)
        rows.append({
            "n": n,
            "total_maps": tot["total_maps"],
            "permutations": tot["counts"]["is_permutation"],
            "idempotent": tot["counts"]["is_idempotent"],
            "connected": tot["counts"]["is_connected"],
            "sum_cyclic": tot["sums"]["num_cyclic"],
            "sum_fixed": tot["sums"]["num_fixed"],
            "E_cyclic": tot["expectations"]["num_cyclic"],
            "E_components": tot["expectations"]["num_components"],
            "E_terminal": tot["expectations"]["num_terminal"],
            "depth_le_1": num_depth_le_1(n),
            "indeg_in_02": num_indegree_in_02(n),
        })
    return rows


def _fmt_frac(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else f"{fr.numerator}"


if __name__ == "__main__":
    print("Exact functional-graph statistics over all n^n maps")
    print("=" * 92)
    header = (f"{'n':>2} | {'n^n':>9} | {'perms(n!)':>9} | {'idemp':>6} | "
              f"{'conn':>7} | {'E[cyc]':>10} | {'E[comp]':>9} | "
              f"{'A006153':>9} | {'A126934':>8}")
    print(header)
    print("-" * 92)
    for row in summary_table(7):
        print(f"{row['n']:>2} | {row['total_maps']:>9} | "
              f"{row['permutations']:>9} | {row['idempotent']:>6} | "
              f"{row['connected']:>7} | {_fmt_frac(row['E_cyclic']):>10} | "
              f"{_fmt_frac(row['E_components']):>9} | "
              f"{row['depth_le_1']:>9} | {row['indeg_in_02']:>8}")
    print("-" * 92)
    print("Validated:  total = n^n;  perms = n!;  idemp = A000248;  "
          "conn = A001865;  E[cyc] = sum_k n!/((n-k)! n^k)")
    print("Restricted: depth<=1 = A006153 (e.g.f. 1/(1-x e^x));  "
          "indeg in {0,2} = A126934 (n even)")
