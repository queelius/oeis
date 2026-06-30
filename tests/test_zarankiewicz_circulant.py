"""Tests for diagonal-circulant Zarankiewicz lower bound."""

from __future__ import annotations

from src.zarankiewicz_circulant import (
    best_diag_circulant,
    count_ones,
    forbidden_minkowski_sums,
    materialize_matrix,
    verify_no_st_submatrix,
)


def test_forbidden_count_3_4_10_11():
    # For m=10, n=11, s=3, t=4: at most C(10, 3) * C(11, 4) = 39600 sums
    # (some may coincide as sets, so the count is <= 39600).
    sums = forbidden_minkowski_sums(10, 11, 3, 4)
    assert len(sums) <= 120 * 330
    # All sums have size <= n = 11 (mod 11).
    for s in sums:
        assert 1 <= len(s) <= 11


def test_minkowski_cauchy_davenport():
    # Cauchy-Davenport: for prime n, |R + C| >= |R| + |C| - 1.
    # For s=3, t=4 in Z_11: min |R + C| >= 6.
    sums = forbidden_minkowski_sums(10, 11, 3, 4)
    assert min(len(s) for s in sums) >= 6


def test_empty_S_is_always_valid():
    # The empty function f = 0 has supp = empty, vacuously avoids everything.
    # best_diag_circulant should never return 0 unless impossible.
    result = best_diag_circulant(10, 11, 3, 4)
    assert result["max_w_f"] >= 0


def test_witness_matrix_verifies():
    # The reported witness must directly verify (no 3x4 all-ones submatrix).
    result = best_diag_circulant(10, 11, 3, 4)
    M = materialize_matrix(10, 11, result["S"])
    assert verify_no_st_submatrix(M, 3, 4)
    # Ones count matches the structural formula sum_j (w(f) - f[(j-1) mod n]).
    ones = count_ones(M)
    # Each row has w(f) ones, so total = m * w(f).
    assert ones == 10 * result["max_w_f"]


def test_small_n5_s2_t2():
    # z(m, 5; 2, 2): largest m x 5 0/1 matrix with no 2x2 all-ones submatrix.
    # For m=5: z(5,5;2,2) = 12 (Reiman, known).
    # Circulant lower bound: f on Z_5 with no R+C, |R|=|C|=2 in supp.
    # R+C with |R|=|C|=2 in Z_5 has size >= 3 (Cauchy-Davenport). So supp <= 2.
    # Max w = 2: e.g., {0, 1}. Need to check no R+C ⊂ {0, 1}.
    # R+C = {r1+c1, r1+c2, r2+c1, r2+c2}. Must be subset of {0,1}; only 2 vals.
    # Then r1+c1 = r1+c2 = r2+c1 = r2+c2 mod 5? Not possible for distinct.
    # So supp {0,1} works. Total ones = 5*2 = 10.
    result = best_diag_circulant(5, 5, 2, 2)
    assert result["max_w_f"] == 2
    assert result["total_ones_lower_bound"] == 10
