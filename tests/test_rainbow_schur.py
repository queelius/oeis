"""
Tests for rainbow (anti-Ramsey) Schur numbers in finite abelian groups.

Prior-art validation (the definition/convention is LOCKED by reproducing
published values):

  Cyclic case -- Bevilacqua, King, Kritschgau, Tait, Tebon, Young,
  "Rainbow numbers for x1 + x2 = k x3 in Z_n", arXiv:1809.04576:
    * rb(Z_2,1) = rb(Z_3,1) = 3            (Remark 1)
    * rb(Z_p,1) = 4 for primes p >= 5      (Theorem 1)
    * rb(Z_n,1) = 2 + sum_i a_i (rb(Z_{p_i},1) - 2)  (Theorem 2)
  The coloring is c : Z_n -> [r] (ALL of Z_n, INCLUDING 0); a rainbow
  triple needs three distinct elements with three distinct colors.

  Grid case -- Fallon, Manhart, Miller, Rehm, Warnberg, Zinnel,
  "Rainbow numbers of [m] x [n] for x1 + x2 = x3", arXiv:2301.10349:
    * rb([m] x [n]) = m + n + 1 for m, n >= 2     (main theorem)

New computations / observations validated here:
  * SAT solver agrees with brute force on every small case.
  * rb(G) depends ONLY on |G| across all abelian groups (order-invariance):
    rb(G) = rb(Z_{|G|}, 1) for every abelian group G (tested through |G|=16).
  * Closed form rb(Z_n,1) = 2 + Omega(n) + Omega_{>=5}(n).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from conftest import requires_sat

from src.rainbow_schur import (
    rb_bruteforce,
    rb_sat,
    schur_triples,
    is_rainbow_free,
    abelian_groups,
    group_rank,
    group_exponent,
    is_cyclic,
)


# ---------------------------------------------------------------------------
# Published closed form for the cyclic case (the validation target).
# ---------------------------------------------------------------------------

def _rb_cyclic_published(n):
    """rb(Z_n, 1) from Theorem 2 of arXiv:1809.04576."""
    from sympy import factorint
    if n == 1:
        return 2  # trivial group: ground set {0}, no solutions, rb=|S|+1
    f = factorint(n)
    rb_prime = lambda p: 3 if p in (2, 3) else 4
    return 2 + sum(a * (rb_prime(p) - 2) for p, a in f.items())


# ---------------------------------------------------------------------------
# VALIDATION GATE 1: reproduce published cyclic rb(Z_n, 1) by brute force.
# ---------------------------------------------------------------------------

class TestCyclicValidationGate:
    """Brute force must reproduce the published cyclic values; this locks
    the definition (color ALL of Z_n, rainbow = 3 distinct elts & colors)."""

    def test_remark1_small(self):
        assert rb_bruteforce([2], include_zero=True) == 3
        assert rb_bruteforce([3], include_zero=True) == 3

    def test_theorem1_primes_ge_5(self):
        for p in (5, 7, 11, 13):
            assert rb_bruteforce([p], include_zero=True) == 4, p

    @pytest.mark.parametrize("n", list(range(2, 12)))
    def test_matches_published_formula(self, n):
        assert rb_bruteforce([n], include_zero=True) == _rb_cyclic_published(n)


# ---------------------------------------------------------------------------
# VALIDATION GATE 2: reproduce the grid formula rb([m]x[n]) = m+n+1.
# (The grid is not a group; we build its triples directly and reuse the
# same rainbow-free machinery, confirming the triple logic is correct.)
# ---------------------------------------------------------------------------

def _grid_triples(m, n):
    elems = [(i, j) for i in range(1, m + 1) for j in range(1, n + 1)]
    idx = {e: t for t, e in enumerate(elems)}
    trip = []
    N = len(elems)
    for a in range(N):
        for b in range(a, N):
            s = (elems[a][0] + elems[b][0], elems[a][1] + elems[b][1])
            k = idx.get(s)
            if k is None:
                continue
            if a == b or k == a or k == b:
                continue
            trip.append((a, b, k))
    return elems, trip


def _rb_grid_brute(m, n):
    from src.rainbow_schur import _set_partitions_by_blocks
    elems, trip = _grid_triples(m, n)
    N = len(elems)
    if not trip:
        return N + 1
    best = 1
    for nb, labels in _set_partitions_by_blocks(N):
        if nb <= best:
            continue
        if is_rainbow_free(labels, trip):
            best = nb
    return best + 1


@pytest.mark.parametrize("m,n", [(2, 2), (2, 3), (2, 4), (3, 3), (2, 5)])
def test_grid_formula(m, n):
    assert _rb_grid_brute(m, n) == m + n + 1


# ---------------------------------------------------------------------------
# SAT solver agrees with brute force (so SAT can be trusted on larger G).
# ---------------------------------------------------------------------------

@requires_sat
@pytest.mark.parametrize("orders", [
    [2], [3], [4], [5], [6], [7], [8], [9],
    [2, 2], [2, 3], [3, 3], [2, 4], [2, 2, 2], [3, 4], [2, 6],
])
def test_sat_matches_bruteforce(orders):
    assert rb_sat(orders, include_zero=True) == rb_bruteforce(orders, include_zero=True)


# ---------------------------------------------------------------------------
# MAIN OBSERVATION: order-invariance. rb(G) = rb(Z_|G|, 1) for all abelian G.
# ---------------------------------------------------------------------------

class TestOrderInvariance:
    """Across all abelian groups of a fixed order, rb is constant and equals
    the published cyclic value. This is the central (conjectural) finding."""

    @pytest.mark.parametrize("order", list(range(2, 13)))
    def test_all_groups_of_order_agree(self, order):
        vals = {rb_bruteforce(list(g), include_zero=True)
                for g in abelian_groups(order)}
        assert len(vals) == 1, (order, vals)
        assert vals.pop() == _rb_cyclic_published(order)

    def test_noncyclic_examples(self):
        # genuinely non-cyclic groups vs their cyclic counterpart of same order
        assert rb_bruteforce([2, 2], include_zero=True) == _rb_cyclic_published(4)
        assert rb_bruteforce([2, 2, 2], include_zero=True) == _rb_cyclic_published(8)
        assert rb_bruteforce([2, 4], include_zero=True) == _rb_cyclic_published(8)
        assert rb_bruteforce([3, 3], include_zero=True) == _rb_cyclic_published(9)

    @requires_sat
    def test_elementary_abelian_order16(self):
        # Z2^4 is maximally far from cyclic Z16 yet gives the same rb.
        assert rb_sat([2, 2, 2, 2], include_zero=True) == _rb_cyclic_published(16)
        assert rb_sat([4, 4], include_zero=True) == _rb_cyclic_published(16)


# ---------------------------------------------------------------------------
# Closed form for the cyclic sequence.
# ---------------------------------------------------------------------------

def test_closed_form_omega():
    """rb(Z_n,1) = 2 + Omega(n) + Omega_{>=5}(n)."""
    from sympy import factorint
    for n in range(2, 500):
        f = factorint(n)
        bigO = sum(f.values())
        o5 = sum(a for p, a in f.items() if p >= 5)
        assert _rb_cyclic_published(n) == 2 + bigO + o5, n


# ---------------------------------------------------------------------------
# Group-invariant helpers.
# ---------------------------------------------------------------------------

def test_group_helpers():
    assert sorted(abelian_groups(16)) == sorted(
        [(16,), (2, 8), (4, 4), (2, 2, 4), (2, 2, 2, 2)])
    assert group_rank((2, 2, 2)) == 3
    assert group_rank((16,)) == 1
    # Z6 in primary-decomposition form is (2, 3): rank 1 over each prime -> cyclic.
    assert is_cyclic((2, 3)) is True
    assert is_cyclic((2, 2)) is False
    assert group_exponent((2, 4)) == 4
    assert group_exponent((3, 3)) == 3
    assert group_exponent((2, 3)) == 6        # Z2 x Z3 ~= Z6
