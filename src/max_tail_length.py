"""
Exact distribution of the MAX TAIL LENGTH (rho-height) of functional graphs.

For an endofunction f : [n] -> [n] the *tail length* of a vertex x is the number
of steps before x enters its cycle (0 for cyclic / recurrent vertices).  The
*max tail length* (a.k.a. rho-height, or simply the *height* of the functional
graph) is

        T(f) = max_{x in [n]} tail_length(x).

This module computes, EXACTLY, the integer count of maps on [n] with a given
max tail length, for every n and t, together with the exact rational mean E[T].

Two equivalent computational routes are provided and cross-checked against each
other and against brute force:

  1. Brute force over all n^n maps (``brute_force_distribution``), reusing the
     per-map analysis in ``finite_map_stats``.  Exact, feasible to n ~ 8.

  2. Exponential generating functions (``height_le_t_counts``), which reach far
     beyond brute force and yield the closed form below.

------------------------------------------------------------------------------
CLOSED FORM (proved; see discoveries/number-theory/max-tail-length/README.md)
------------------------------------------------------------------------------
A functional graph is a SET of CYCLES of ROOTED TREES (Flajolet-Sedgewick,
Analytic Combinatorics II.4): the recurrent vertices form the cyclic core (any
permutation of them), and a rooted tree hangs off each cyclic vertex.  A vertex
at depth d inside such a tree (root = cyclic vertex at depth 0) has tail length
exactly d.  Hence

        T(f) <= t   <=>   every tree attached to the core has height <= t,

where a rooted tree has *height* = max depth of a vertex (a lone root: 0).

Let e_t(x) be the EGF for rooted trees of height <= t.  Removing the root and
recursing on the children-forest gives the recursion

        e_0(x) = x,        e_t(x) = x * exp(e_{t-1}(x))        (t >= 1).

A functional graph whose every tree has height <= t is a set of cycles of such
trees, with EGF (SET o CYCLE o (trees of height <= t)):

        F_{<=t}(x) = 1 / (1 - e_t(x)).

Therefore the number of maps on [n] with T <= t is

        H(n, t) = n! * [x^n]  1 / (1 - e_t(x)),                       (CUMULATIVE)

and the number with T = t exactly is the column difference

        D(n, t) = H(n, t) - H(n, t-1)
                = n! * [x^n] ( 1/(1 - e_t) - 1/(1 - e_{t-1}) ).        (EXACT-t)

Boundary / sanity identities (all proved, all tested):
  * e_{-1} := 0, so H(n,-1) = [n = 0] and D(n,0) = H(n,0) = n!  (permutations:
    T = 0 iff every vertex is cyclic iff f is a bijection).
  * H(n, t) = n^n for all t >= n-1  (the longest possible tail on [n] is n-1).
  * D(n, n-1) = n!  (the height-(n-1) maps: a single path of n vertices whose
    end is a fixed point -- n! labelings).
  * Row sums  sum_t D(n,t) = n^n.

------------------------------------------------------------------------------
OEIS
------------------------------------------------------------------------------
  * D(n,t) read by rows is OEIS A216242 (Geoffrey Critzer, Mar 14 2013),
    "number of functions f:{1,...,n}->{1,...,n} with a height of k".  Critzer's
    published e.g.f. for column k, 1/(1-G(k+1)) - 1/(1-G(k)) with
    G(k) = x*exp(G(k-1)), G(0) = 0, is exactly the formula above under the
    indexing G(k) = e_{k-1}.  This module reproduces every term of A216242 and
    the present derivation supplies a self-contained structural proof.
  * The t = 1 cumulative column H(n,1) = 1, 4, 21, 148, 1305, ... is OEIS
    A006153 (e.g.f. 1/(1 - x e^x)), the "distance <= 1 from a cycle" maps.
  * The full CUMULATIVE triangle H(n,t) is (as of 2026-06) NOT in OEIS.

------------------------------------------------------------------------------
ASYMPTOTICS
------------------------------------------------------------------------------
Flajolet & Odlyzko (1990) show the expected height of a random mapping on [n]
satisfies  E[T] ~ c * sqrt(n)  with  c = sqrt(pi/2)*... ; concretely the mean
is asymptotic to a constant multiple of sqrt(n) of the same order as
sqrt(pi n / 2).  ``expected_max_tail`` returns the exact Fraction and the trend
is checked numerically against sqrt(pi n / 2) in the tests.

References:
  - Flajolet & Odlyzko (1990), "Random Mapping Statistics", EUROCRYPT '89,
    LNCS 434, pp. 329-354 (height of random mappings).
  - Flajolet & Sedgewick (2009), "Analytic Combinatorics", Cambridge UP, II.4.
  - OEIS A216242 (Critzer), A006153, A000312 (n^n), A000142 (n!).
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction
from typing import Iterator

from src.finite_map_stats import all_maps, tail_lengths


# ---------------------------------------------------------------------------
# Per-map max tail length and brute-force distribution
# ---------------------------------------------------------------------------

def max_tail_length(f: tuple[int, ...] | list[int]) -> int:
    """Max tail length (rho-height) of a single endofunction f on [n].

    0 iff f is a permutation (every vertex cyclic); the maximum over [n] of the
    number of steps each vertex takes to reach its cycle.
    """
    tails = tail_lengths(f)
    return max(tails) if tails else 0


def _max_tail_fast(f: tuple[int, ...], n: int) -> int:
    """Max tail length of f on [n], optimized for the brute-force inner loop.

    Single-pass memoized depth computation.  state[x]: 0 unseen, 1 on the path
    currently being explored, 2 resolved (depth known in ``depth``).  Cyclic
    vertices get depth 0; a tail vertex gets 1 + depth of its image.  Returns
    the maximum depth.  Equivalent to ``max_tail_length`` but ~10x faster.
    """
    depth = [-1] * n          # resolved tail length, -1 if unknown
    state = [0] * n           # 0 unseen, 1 in-progress, 2 done
    best = 0
    for s in range(n):
        if state[s]:
            continue
        path = []
        x = s
        while state[x] == 0:
            state[x] = 1
            path.append(x)
            x = f[x]
        if state[x] == 1:
            # hit a vertex on the current path => found a fresh cycle; every
            # vertex from x onward in ``path`` is cyclic (depth 0).
            ci = path.index(x)
            for y in path[ci:]:
                depth[y] = 0
                state[y] = 2
            base = 0            # x is cyclic, depth 0
            tail = path[:ci]
        else:
            # x already resolved (state 2): its depth is known.
            base = depth[x]
            tail = path
        # assign depths to the tail vertices, walking back from x.
        d = base
        for y in reversed(tail):
            d += 1
            depth[y] = d
            state[y] = 2
            if d > best:
                best = d
    return best


def brute_force_distribution(n: int) -> dict[int, int]:
    """Exact distribution {t: count} of the max tail length over all n^n maps.

    Pure exhaustion -- exact integers, the ground truth used to validate the
    EGF route.  Feasible up to n = 8 (8^8 = 16,777,216 maps).
    """
    dist: Counter[int] = Counter()
    for f in all_maps(n):
        dist[_max_tail_fast(f, n)] += 1
    return dict(dist)


# ---------------------------------------------------------------------------
# Exact truncated EGF arithmetic over the rationals
# ---------------------------------------------------------------------------
# A series is represented by its coefficient list [a_0, a_1, ..., a_N] of the
# ORDINARY coefficients a_k = [x^k] (NOT divided by k!).  All arithmetic is
# exact (Fraction).  We only need: multiply, the geometric series 1/(1-g) for a
# series g with g(0)=0, and exp(g) for a series g with g(0)=0.

def _mul(a: list[Fraction], b: list[Fraction], N: int) -> list[Fraction]:
    """Truncated product of two power series, keeping degrees 0..N."""
    out = [Fraction(0)] * (N + 1)
    for i, ai in enumerate(a):
        if ai == 0 or i > N:
            continue
        for j, bj in enumerate(b):
            if j + i > N:
                break
            if bj:
                out[i + j] += ai * bj
    return out


def _exp_zero_const(g: list[Fraction], N: int) -> list[Fraction]:
    """exp(g) for a power series g with g(0) = 0, truncated to degree N.

    Uses the ODE  E' = g' E  =>  n e_n = sum_{k=1..n} k g_k e_{n-k}, e_0 = 1.
    Exact rational arithmetic throughout.
    """
    assert not g or g[0] == 0, "exp requires zero constant term"
    e = [Fraction(0)] * (N + 1)
    e[0] = Fraction(1)
    for m in range(1, N + 1):
        s = sum(k * (g[k] if k < len(g) else Fraction(0)) * e[m - k]
                for k in range(1, m + 1))
        e[m] = s / m
    return e


def _one_over_one_minus(g: list[Fraction], N: int) -> list[Fraction]:
    """1 / (1 - g) for a power series g with g(0) = 0, truncated to degree N.

    Geometric series via the recurrence h = 1 + g*h, i.e.
    h_0 = 1, h_m = sum_{k=1..m} g_k h_{m-k}.
    """
    assert not g or g[0] == 0, "1/(1-g) needs g(0)=0 for a power series"
    h = [Fraction(0)] * (N + 1)
    h[0] = Fraction(1)
    for m in range(1, N + 1):
        h[m] = sum((g[k] if k < len(g) else Fraction(0)) * h[m - k]
                   for k in range(1, m + 1))
    return h


def tree_egf_height_le(t: int, N: int) -> list[Fraction]:
    """EGF coefficients (ordinary [x^k], not /k!) of rooted trees of height <= t,
    truncated to degree N.

    e_0 = x;  e_t = x * exp(e_{t-1}).  Convention: t = -1 returns the zero series
    (no tree), which makes the cumulative formula H(n,-1) = [n=0] come out right.
    """
    if t < 0:
        return [Fraction(0)] * (N + 1)
    e = [Fraction(0)] * (N + 1)
    if N >= 1:
        e[1] = Fraction(1)  # e_0 = x
    x = [Fraction(0)] * (N + 1)
    if N >= 1:
        x[1] = Fraction(1)
    for _ in range(t):
        e = _mul(x, _exp_zero_const(e, N), N)
    return e


def _egf_to_counts(series: list[Fraction], N: int) -> list[int]:
    """Convert an EGF (ordinary coefficients a_k) to labeled counts n! a_n."""
    out = []
    for n in range(N + 1):
        c = series[n] * math.factorial(n)
        assert c.denominator == 1, f"non-integral count at n={n}: {c}"
        out.append(int(c))
    return out


# ---------------------------------------------------------------------------
# Cumulative and exact-t counts via EGF (the closed form)
# ---------------------------------------------------------------------------

def height_le_t_counts(t: int, N: int) -> list[int]:
    """H(n, t) = number of maps on [n] with max tail length <= t, for n=0..N.

    Closed form:  H(n, t) = n! [x^n] 1/(1 - e_t(x)).  Returns the list indexed
    by n.  For t < 0 this is [1, 0, 0, ...] (only the empty map has height <= -1).
    """
    e = tree_egf_height_le(t, N)
    F = _one_over_one_minus(e, N)
    return _egf_to_counts(F, N)


def height_eq_t_counts(t: int, N: int) -> list[int]:
    """D(n, t) = number of maps on [n] with max tail length EXACTLY t, n=0..N.

    D(n, t) = H(n, t) - H(n, t-1).  Equivalently the EGF column difference
    1/(1-e_t) - 1/(1-e_{t-1}).
    """
    upper = height_le_t_counts(t, N)
    lower = height_le_t_counts(t - 1, N)
    return [upper[n] - lower[n] for n in range(N + 1)]


def distribution_via_egf(n: int) -> dict[int, int]:
    """Exact {t: count} distribution of max tail length on [n], via the EGF.

    Only t in 0..n-1 can occur for n >= 1 (the empty map n=0 has height 0).
    """
    if n == 0:
        return {0: 1}
    dist: dict[int, int] = {}
    cum_prev = 0
    H = [height_le_t_counts(t, n)[n] for t in range(n)]  # H(n,0..n-1)
    for t in range(n):
        d = H[t] - cum_prev
        if d:
            dist[t] = d
        cum_prev = H[t]
    return dist


def distribution_triangle(n_max: int) -> list[list[int]]:
    """Rows D(n, 0..n-1) for n = 1..n_max -- the OEIS A216242 triangle."""
    rows = []
    # Compute all columns up to n_max-1 once, sliced per n.
    cum = {t: height_le_t_counts(t, n_max) for t in range(-1, n_max)}
    for n in range(1, n_max + 1):
        row = [cum[t][n] - cum[t - 1][n] for t in range(n)]
        rows.append(row)
    return rows


def cumulative_triangle(n_max: int) -> list[list[int]]:
    """Rows H(n, 0..n-1) for n = 1..n_max -- cumulative counts (NOT in OEIS)."""
    cum = {t: height_le_t_counts(t, n_max) for t in range(n_max)}
    rows = []
    for n in range(1, n_max + 1):
        rows.append([cum[t][n] for t in range(n)])
    return rows


# ---------------------------------------------------------------------------
# Expected max tail length (exact rational) and closed-form helpers
# ---------------------------------------------------------------------------

def expected_max_tail(n: int) -> Fraction:
    """Exact mean E[T] of the max tail length over a uniform random map on [n].

    E[T] = (1/n^n) sum_t t * D(n, t).  A numerically convenient equivalent uses
    the tail-sum identity  E[T] = (1/n^n) sum_{t>=1} (n^n - H(n, t-1)).
    Returned as an exact Fraction.
    """
    if n == 0:
        return Fraction(0)
    nn = n ** n
    # E[T] = sum_{t>=1} P(T >= t) = sum_{t>=1} (n^n - H(n,t-1))/n^n,
    # summed over t=1..n-1 (P(T>=t)=0 for t>=n).
    total = 0
    for t in range(1, n):
        total += nn - height_le_t_counts(t - 1, n)[n]
    return Fraction(total, nn)


def num_permutations(n: int) -> int:
    """n! -- equals D(n, 0) = H(n, 0), the maps with max tail length 0."""
    return math.factorial(n)


def closed_form_diagonal_top(n: int) -> int:
    """D(n, n-1) = n! -- closed form for the largest possible height on [n]."""
    return math.factorial(n)


# ---------------------------------------------------------------------------
# Summary driver
# ---------------------------------------------------------------------------

def _fmt_frac(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else f"{fr.numerator}"


def summary(n_max: int = 10) -> None:
    """Print the exact-T distribution triangle, sanity checks, and E[T] trend."""
    print("Max tail length (rho-height) of functional graphs -- exact counts")
    print("=" * 78)
    print("Exact-T distribution triangle D(n,t), t = 0..n-1  (OEIS A216242):")
    for n, row in enumerate(distribution_triangle(n_max), start=1):
        rs = sum(row)
        flag = "OK" if rs == (n ** n) and row[0] == math.factorial(n) \
            and row[-1] == math.factorial(n) else "??"
        print(f"  n={n:2d}: {row}   (sum={rs}=n^n {flag})")
    print()
    print("E[T] trend vs sqrt(pi n / 2):")
    print(f"  {'n':>3} | {'E[T] exact':>22} | {'E[T]':>9} | {'sqrt(pi n/2)':>12}")
    for n in range(1, n_max + 1):
        et = expected_max_tail(n)
        print(f"  {n:>3} | {_fmt_frac(et):>22} | {float(et):>9.5f} | "
              f"{math.sqrt(math.pi * n / 2):>12.5f}")


if __name__ == "__main__":
    summary(10)
