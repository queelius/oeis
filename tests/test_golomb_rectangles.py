"""Tests for the Golomb-rectangle solver G(m,n).

Layers:
1. Geometry / property predicate is correct (brute checks).
2. The dihedral symmetry group is exactly the set of grid motions that
   preserve the distinct-difference property (and double-lex is NOT, which is
   why we use dihedral lex-leader symmetry breaking, not double-lex).
3. SOUNDNESS: SAT (with and without symmetry breaking) reproduces the
   brute-force ground truth, cell by cell, for all small m x n.
4. Prior-art cross-check: the m=1 boundary row equals OEIS A143824
   (1D Golomb-ruler / Sidon-set sizes), and G(m,n) is symmetric.
"""

from __future__ import annotations

import itertools

import pytest

from src.golomb_rectangles import (
    KNOWN_G,
    brute_force_max,
    cells,
    dihedral_images,
    golomb_rectangle,
    sat_can_place,
    verify_witness,
    _bad_quadruples,
    _has_distinct_differences,
)

pytest.importorskip("pysat")


# Brute-force ground truth used across tests (cheap to recompute).
BRUTE = {(m, n): brute_force_max(m, n)[0]
         for m in range(1, 5) for n in range(m, 6)}


# ---- (1) property predicate ----

def test_property_predicate_basic():
    # A single dot and the empty set are trivially Golomb.
    assert _has_distinct_differences([])
    assert _has_distinct_differences([(0, 0)])
    # 3 dots, all differences distinct.
    assert _has_distinct_differences([(0, 0), (0, 1), (1, 0)])
    # Square's 4 corners: (1,1)-(0,0)=(1,1) and (1,0)-(0,-1)... actually the
    # 2x2 full grid repeats the difference (1,0) and (0,1), so NOT Golomb.
    assert not _has_distinct_differences([(0, 0), (0, 1), (1, 0), (1, 1)])


def test_known_small_values_match_brute():
    for (m, n), g in KNOWN_G.items():
        assert brute_force_max(m, n)[0] == g, (m, n)


@pytest.mark.parametrize("m,n,rmax", [(3, 3, 9), (3, 4, 8), (2, 5, 8)])
def test_bad_quadruples_exactly_characterize_property(m, n, rmax):
    # CORE soundness: a set of cells violates the distinct-difference property
    # IFF it contains some forbidden "bad quadruple" as a subset.  This is the
    # exact statement the CNF clauses encode.
    bad = [set(q) for q in _bad_quadruples(m, n)]
    grid = cells(m, n)
    # Forbidden sets only have size 3 or 4 (size 2 impossible: needs 2d=0).
    assert set(len(b) for b in bad) <= {3, 4}
    for r in range(2, rmax + 1):
        for S in itertools.combinations(grid, r):
            Sset = set(S)
            violates = not _has_distinct_differences(S)
            covered = any(b <= Sset for b in bad)
            assert violates == covered, (m, n, S)


# ---- (2) symmetry group correctness ----

@pytest.mark.parametrize("m,n", [(2, 3), (3, 3), (4, 4), (3, 5), (4, 5)])
def test_dihedral_images_preserve_property_and_stay_on_grid(m, n):
    g, w = brute_force_max(m, n)
    for name, mp in dihedral_images(m, n):
        img = [mp(r, c) for (r, c) in w]
        assert all(0 <= r < m and 0 <= c < n for (r, c) in img), name
        assert _has_distinct_differences(img), \
            f"dihedral image {name} broke the Golomb property"


def test_square_has_eight_symmetries_rect_has_four():
    assert len(dihedral_images(4, 4)) == 8
    assert len(dihedral_images(3, 5)) == 4


def test_nondihedral_row_swap_can_break_property():
    # The crucial fact that rules out double-lex: a plain row permutation need
    # NOT preserve the distinct-difference property.
    g, w = brute_force_max(4, 4)
    assert g == 6
    perm = [1, 0, 2, 3]  # swap rows 0,1 only -- not a dihedral motion
    w2 = [(perm[r], c) for (r, c) in w]
    assert not _has_distinct_differences(w2)


# ---- (3) SAT soundness vs brute force ----

@pytest.mark.parametrize("m,n", sorted(BRUTE))
def test_sat_decision_oracle_matches_brute_all_k(m, n):
    g = BRUTE[(m, n)]
    for k in range(1, m * n + 1):
        expect = (k <= g)
        s_sym, wsym = sat_can_place(m, n, k, symmetry=True)
        s_no, _ = sat_can_place(m, n, k, symmetry=False)
        assert s_sym == expect, (m, n, k, "sym")
        assert s_no == expect, (m, n, k, "nosym")
        if s_sym:
            assert verify_witness(m, n, wsym)


@pytest.mark.parametrize("m,n", sorted(BRUTE))
def test_golomb_rectangle_value_and_proof(m, n):
    res_sym = golomb_rectangle(m, n, symmetry=True)
    res_no = golomb_rectangle(m, n, symmetry=False)
    assert res_sym["G"] == res_no["G"] == BRUTE[(m, n)], (m, n)
    assert res_sym["proven"] and res_no["proven"]
    assert verify_witness(m, n, res_sym["witness"])
    assert verify_witness(m, n, res_no["witness"])


# ---- (4) prior-art cross-check ----

# OEIS A143824(n): size of largest subset of {1..n} with all distinct
# differences (a(0)=0).  G(1,n) (single row of length n) must equal A143824(n).
A143824 = [0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5]


@pytest.mark.parametrize("n", list(range(1, 13)))
def test_single_row_equals_A143824(n):
    assert brute_force_max(1, n)[0] == A143824[n], n


@pytest.mark.parametrize("m,n", [(2, 3), (2, 5), (3, 4), (3, 5), (4, 6)])
def test_symmetry_of_G(m, n):
    assert golomb_rectangle(m, n)["G"] == golomb_rectangle(n, m)["G"]
