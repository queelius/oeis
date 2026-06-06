"""
BIVARIATE distribution of functional graphs refining the max-tail-length
(rho-height) by a second statistic: either the number of CYCLIC POINTS or the
number of COMPONENTS.

Background (univariate, see src/max_tail_length.py).  For an endofunction
f : [n] -> [n] the *max tail length* T(f) is the longest distance from any
vertex to its cycle (0 iff f is a permutation).  The univariate count
D(n,t) = #{f : T(f) = t} read by rows is OEIS A216242 (Critzer 2013), with EGF
D(n,t) = n! [x^n] (1/(1-e_t) - 1/(1-e_{t-1})), where e_t is the EGF of rooted
trees of height <= t:  e_0 = x,  e_t = x * exp(e_{t-1}).

This module computes, EXACTLY, the two JOINT distributions

    B(n,t,c) = #{ f : T(f) = t  AND  f has exactly c cyclic points },
    Bm(n,t,m) = #{ f : T(f) = t  AND  f has exactly m components },

for all n and t and (c or m), validates them against brute force at small n,
validates both marginals against the literature, and supplies a proved
bivariate EGF for each.

------------------------------------------------------------------------------
BIVARIATE EGFs (proved; see discoveries/.../bivariate-README.md)
------------------------------------------------------------------------------
A functional graph is a SET of CYCLES of ROOTED TREES (Flajolet-Sedgewick
II.4): the cyclic points carry a permutation (the core), and a rooted tree
hangs off each cyclic point with the cyclic point as its root.  A vertex at
depth d in such a tree has tail length d, so

        T(f) <= t   <=>   every hanging tree has height <= t.                (★)

Let e_t(x) be the EGF of rooted trees of height <= t.

* CYCLIC POINTS, marked by u.  The core is a permutation, i.e. a SEQUENCE of
  "slots"; each slot is filled by one rooted tree of height <= t whose root is
  a cyclic point.  Marking that cyclic point by u replaces each slot's weight
  e_t by u*e_t, and a sequence of slots has EGF 1/(1 - (u e_t)).  Hence

        F_{<=t}(x,u) = 1 / (1 - u * e_t(x)),                            (CUM-c)
        B(n,t,c)     = n! [x^n u^c] ( 1/(1 - u e_t) - 1/(1 - u e_{t-1}) ).

* COMPONENTS, marked by u.  A component is one cycle of trees, EGF
  L_t = log(1/(1 - e_t)).  A functional graph is a SET of components, so
  marking each component by u gives exp(u * L_t):

        G_{<=t}(x,u) = exp( u * log(1/(1 - e_t(x))) ) = (1 - e_t(x))^{-u},
        Bm(n,t,m)    = n! [x^n u^m] ( G_{<=t} - G_{<=t-1} ).            (CUM-m)

------------------------------------------------------------------------------
MARGINALS (both validated exactly against the literature)
------------------------------------------------------------------------------
  * sum_c B(n,t,c) = D(n,t)               (A216242; univariate max tail length)
  * sum_t B(n,t,c) = M(n,c) = C(n,c) c! * c * n^{n-c-1}   (c<n; = 1 forest for
        c=n), the classical "endofunctions with c cyclic points", OEIS A066324.
  * sum_m Bm(n,t,m) = D(n,t)              (A216242)
  * sum_t Bm(n,t,m) = A060281(n,m)        (endofunctions with m components/cycles)

------------------------------------------------------------------------------
OEIS STATUS (searched Jun 2026; OEIS query via curl, WebFetch 403-blocked)
------------------------------------------------------------------------------
  * Both MARGINAL triangles are classical: A066324 (cyclic points) and A060281
    (components); the univariate max-tail marginal is A216242.
  * Neither JOINT B(n,t,c) nor Bm(n,t,m) was found in OEIS in any natural row
    reading (full-square, support-only, or by (n,c,t)).  They appear to be new;
    a submission draft for the (t,c) joint lives in the discoveries directory.

References:
  - P. Flajolet, A. Odlyzko, "Random Mapping Statistics", EUROCRYPT '89,
    LNCS 434 (1990), 329-354.
  - P. Flajolet, R. Sedgewick, "Analytic Combinatorics", Cambridge UP (2009),
    II.4 (functional graphs as sets of cycles of trees), III.7 (markings).
  - OEIS A216242 (Critzer; max tail length), A066324 (cyclic points),
    A060281 (components / cycles), A000312 (n^n).
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction
from itertools import product

from src.finite_map_stats import cyclic_points, num_components
from src.max_tail_length import _max_tail_fast, tree_egf_height_le


# ---------------------------------------------------------------------------
# Brute-force joint distributions (ground truth)
# ---------------------------------------------------------------------------

def brute_force_joint_tc(n: int) -> dict[tuple[int, int], int]:
    """Exact joint distribution {(t, c): count} over all n^n maps, where t is
    the max tail length and c the number of cyclic points.  Feasible to n ~ 8.
    """
    j: Counter[tuple[int, int]] = Counter()
    for f in product(range(n), repeat=n):
        j[(_max_tail_fast(f, n), len(cyclic_points(f)))] += 1
    return dict(j)


def brute_force_joint_tm(n: int) -> dict[tuple[int, int], int]:
    """Exact joint distribution {(t, m): count} over all n^n maps, where t is
    the max tail length and m the number of (weakly-connected) components."""
    j: Counter[tuple[int, int]] = Counter()
    for f in product(range(n), repeat=n):
        j[(_max_tail_fast(f, n), num_components(f))] += 1
    return dict(j)


# ---------------------------------------------------------------------------
# Exact bivariate truncated power-series arithmetic
# ---------------------------------------------------------------------------
# A bivariate series sum_{i,j} a_{i,j} x^i u^j is stored as a list indexed by
# the x-degree i (0..N); entry i is a dict {j: Fraction} of the u-coefficients.
# x is the labelled variable (we divide by i! only at the very end); u is an
# ordinary marking variable (its powers are NOT factorial-weighted).

def _bz(N: int) -> list[dict[int, Fraction]]:
    """The zero bivariate series, x-degrees 0..N."""
    return [dict() for _ in range(N + 1)]


def _badd_inplace(dst: dict[int, Fraction], j: int, val: Fraction) -> None:
    if val:
        dst[j] = dst.get(j, Fraction(0)) + val


def _univariate_to_bivariate(e: list[Fraction], u_power: int,
                             N: int) -> list[dict[int, Fraction]]:
    """Lift a plain x-series e (list of Fraction) to a bivariate series whose
    every term carries u^{u_power}: x^i coeff e_i -> e_i x^i u^{u_power}."""
    out = _bz(N)
    for i in range(min(len(e), N + 1)):
        if e[i]:
            out[i] = {u_power: e[i]}
    return out


def _geom_in_ue(ue: list[dict[int, Fraction]], N: int) -> list[dict[int, Fraction]]:
    """1 / (1 - g) for a bivariate series g with g(0,*) = 0 (no constant term),
    via h = 1 + g*h:  h_0 = 1 (at u^0),  h_i = sum_{k=1..i} g_k * h_{i-k}."""
    h = _bz(N)
    h[0] = {0: Fraction(1)}
    for i in range(1, N + 1):
        acc: dict[int, Fraction] = {}
        for k in range(1, i + 1):
            gk = ue[k]
            if not gk:
                continue
            hik = h[i - k]
            for ju, ca in gk.items():
                for ku, cb in hik.items():
                    _badd_inplace(acc, ju + ku, ca * cb)
        h[i] = acc
    return h


def _log_one_over_one_minus(e: list[Fraction], N: int) -> list[Fraction]:
    """L = log(1/(1-e)) = sum_{k>=1} e^k / k for a plain x-series e with
    e(0) = 0, truncated to x-degree N.  Exact Fractions."""
    L = [Fraction(0)] * (N + 1)
    ek = [Fraction(0)] * (N + 1)
    ek[0] = Fraction(1)  # e^0

    def mul(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
        o = [Fraction(0)] * (N + 1)
        for i, pi in enumerate(p):
            if not pi:
                continue
            for jj, qj in enumerate(q):
                if i + jj > N:
                    break
                if qj:
                    o[i + jj] += pi * qj
        return o

    for k in range(1, N + 1):
        ek = mul(ek, e)
        if not any(ek):
            break
        for i in range(N + 1):
            if ek[i]:
                L[i] += ek[i] / k
    return L


def _exp_u_times(L: list[Fraction], N: int) -> list[dict[int, Fraction]]:
    """exp(u * L) for a plain x-series L with L(0) = 0, returned as a bivariate
    series in (x, u).  Since L has x-valuation >= 1, L^m has x-valuation >= m,
    so the sum_{m>=0} u^m L^m / m! truncates at m = N."""
    G = _bz(N)
    G[0] = {0: Fraction(1)}
    Lm = [Fraction(0)] * (N + 1)
    Lm[0] = Fraction(1)  # L^0

    def mul(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
        o = [Fraction(0)] * (N + 1)
        for i, pi in enumerate(p):
            if not pi:
                continue
            for jj, qj in enumerate(q):
                if i + jj > N:
                    break
                if qj:
                    o[i + jj] += pi * qj
        return o

    for m in range(1, N + 1):
        Lm = mul(Lm, L)
        if not any(Lm):
            break
        fac = math.factorial(m)
        for i in range(N + 1):
            if Lm[i]:
                _badd_inplace(G[i], m, Lm[i] / fac)
    return G


# ---------------------------------------------------------------------------
# Cumulative bivariate EGFs and exact-t joint counts
# ---------------------------------------------------------------------------

def cum_tc_egf(t: int, N: int) -> list[dict[int, Fraction]]:
    """Bivariate EGF F_{<=t}(x,u) = 1/(1 - u e_t(x)) as a (x,u) series to
    x-degree N.  Coefficient n! [x^n u^c] is #{f : T(f) <= t, c cyclic points}.
    For t < 0 returns the empty-map series (constant 1)."""
    e = tree_egf_height_le(t, N)            # e_{-1} = 0 -> F = 1
    ue = _univariate_to_bivariate(e, 1, N)  # u * e_t
    return _geom_in_ue(ue, N)


def cum_tm_egf(t: int, N: int) -> list[dict[int, Fraction]]:
    """Bivariate EGF G_{<=t}(x,u) = exp(u * log(1/(1-e_t(x)))) = (1-e_t)^{-u}
    to x-degree N.  Coefficient n! [x^n u^m] is #{f : T(f) <= t, m components}."""
    e = tree_egf_height_le(t, N)
    L = _log_one_over_one_minus(e, N)       # e_{-1}=0 -> L=0 -> G=1
    return _exp_u_times(L, N)


def _coeff_counts(cum_egf, t: int, n: int) -> dict[int, int]:
    """n! [x^n u^j] ( cum(t) - cum(t-1) ) as {j: integer count}, dropping zeros."""
    nfac = math.factorial(n)
    cur = cum_egf(t, n)[n]
    prev = cum_egf(t - 1, n)[n]
    out: dict[int, int] = {}
    for j in set(cur) | set(prev):
        val = (cur.get(j, Fraction(0)) - prev.get(j, Fraction(0))) * nfac
        assert val.denominator == 1, f"non-integral joint count n={n} t={t} j={j}: {val}"
        iv = int(val)
        if iv:
            out[j] = iv
    return out


def joint_tc_via_egf(n: int) -> dict[tuple[int, int], int]:
    """Exact joint {(t, c): count} for [n] via the bivariate EGF 1/(1-u e_t)."""
    out: dict[tuple[int, int], int] = {}
    for t in range(n):
        for c, v in _coeff_counts(cum_tc_egf, t, n).items():
            out[(t, c)] = v
    return out


def joint_tm_via_egf(n: int) -> dict[tuple[int, int], int]:
    """Exact joint {(t, m): count} for [n] via the bivariate EGF (1-e_t)^{-u}."""
    out: dict[tuple[int, int], int] = {}
    for t in range(n):
        for m, v in _coeff_counts(cum_tm_egf, t, n).items():
            out[(t, m)] = v
    return out


# ---------------------------------------------------------------------------
# Marginals (closed forms) for validation
# ---------------------------------------------------------------------------

def cyclic_points_marginal(n: int) -> list[int]:
    """M(n,c) for c = 1..n: number of endofunctions on [n] with exactly c
    cyclic points = C(n,c) * c! * (rooted forests on n nodes, c given roots).

    The forest count is c * n^{n-c-1} for c < n (generalized Cayley) and 1 for
    c = n (the all-roots forest).  Row n of OEIS A066324.
    """
    out = []
    for c in range(1, n + 1):
        forests = 1 if c == n else c * (n ** (n - c - 1))
        out.append(math.comb(n, c) * math.factorial(c) * forests)
    return out


# ---------------------------------------------------------------------------
# Triangles (read by rows) for OEIS / display
# ---------------------------------------------------------------------------

def joint_tc_rows(n: int, *, full_square: bool = False) -> list[list[int]]:
    """Rows of B(n, t, c) for fixed n.  Row index is t = 0..n-1; within a row
    the columns are c.

    full_square=False (default): support-only.  For t = 0 the only nonzero
    entry is c = n (a permutation), emitted as a length-1 row [n!]; for t >= 1
    the columns are c = 1..n-t (a max tail of t needs t non-cyclic vertices).

    full_square=True: every row has columns c = 1..n with explicit zeros.
    """
    j = joint_tc_via_egf(n)
    rows = []
    for t in range(n):
        if full_square:
            rows.append([j.get((t, c), 0) for c in range(1, n + 1)])
        elif t == 0:
            rows.append([j.get((0, n), 0)])
        else:
            rows.append([j.get((t, c), 0) for c in range(1, n - t + 1)])
    return rows


def joint_tm_rows(n: int, *, full_square: bool = False) -> list[list[int]]:
    """Rows of Bm(n, t, m) for fixed n.  Row index t = 0..n-1; columns m.

    Support-only: t = 0 gives m = 1..n (a permutation can have any number of
    cycles); t >= 1 gives m = 1..n-t.  full_square pads to m = 1..n with zeros.
    """
    j = joint_tm_via_egf(n)
    rows = []
    for t in range(n):
        if full_square:
            rows.append([j.get((t, m), 0) for m in range(1, n + 1)])
        elif t == 0:
            rows.append([j.get((0, m), 0) for m in range(1, n + 1)])
        else:
            rows.append([j.get((t, m), 0) for m in range(1, n - t + 1)])
    return rows


def flatten_joint_tc_support(n_max: int) -> list[int]:
    """The (t,c) joint flattened for OEIS: for n = 1..n_max, read rows
    t = 0,1,..,n-1, and within each row the nonzero columns c (the support).
    """
    out: list[int] = []
    for n in range(1, n_max + 1):
        for row in joint_tc_rows(n):
            out.extend(row)
    return out


def flatten_joint_tm_support(n_max: int) -> list[int]:
    """The (t,m) joint flattened for OEIS (support-only, by n then t then m)."""
    out: list[int] = []
    for n in range(1, n_max + 1):
        for row in joint_tm_rows(n):
            out.extend(row)
    return out


# ---------------------------------------------------------------------------
# Summary driver
# ---------------------------------------------------------------------------

def summary(n_max: int = 7) -> None:
    """Print the joint (t,c) and (t,m) tables, marginal validations, and the
    flattened OEIS-candidate sequences."""
    print("Bivariate refinement of the max-tail-length distribution")
    print("=" * 74)
    print("\nJoint B(n,t,c): rows t=0..n-1, columns c=1..n (full square)")
    for n in range(1, n_max + 1):
        print(f"  n={n}:")
        for t, row in enumerate(joint_tc_rows(n, full_square=True)):
            print(f"    t={t}: {row}")
    print("\nMarginal checks (must all be True):")
    from src.max_tail_length import distribution_triangle
    Dtri = distribution_triangle(n_max)
    for n in range(1, n_max + 1):
        jtc = joint_tc_via_egf(n)
        Dt = Counter()
        for (t, c), v in jtc.items():
            Dt[t] += v
        d_ok = [Dt[t] for t in range(n)] == Dtri[n - 1]
        Mc = Counter()
        for (t, c), v in jtc.items():
            Mc[c] += v
        m_ok = [Mc[c] for c in range(1, n + 1)] == cyclic_points_marginal(n)
        print(f"  n={n}: sum_c B = A216242 row: {d_ok};  "
              f"sum_t B = A066324 row: {m_ok}")
    print("\nFlattened (t,c) joint, support-only (OEIS candidate), n=1..%d:" % n_max)
    print("  " + ",".join(map(str, flatten_joint_tc_support(n_max))))
    print("\nFlattened (t,m) joint, support-only, n=1..%d:" % n_max)
    print("  " + ",".join(map(str, flatten_joint_tm_support(n_max))))


if __name__ == "__main__":
    summary(7)
