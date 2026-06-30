"""
Tests for the exact max-tail-length (rho-height) distribution of functional
graphs.

Validation strategy (exact arithmetic throughout):

  1. The EGF route (1/(1-e_t), e_t = x*exp(e_{t-1})) is checked term-by-term
     against pure BRUTE FORCE over all n^n maps for n = 1..7 (and n=8 under the
     slow marker).  This is the core 1:1 ground-truth check.

  2. Sanity / boundary identities, all proved in the module docstring:
       * D(n, 0) = H(n, 0) = n!         (max tail length 0  <=>  permutation)
       * D(n, n-1) = n!                  (the longest height = a single path)
       * row sums  sum_t D(n,t) = n^n
       * H(n, t) = n^n for all t >= n-1

  3. Prior-art / OEIS validation:
       * D(n,t) read by rows == OEIS A216242 (Critzer 2013), published DATA.
       * H(n,1) (height <= 1) == OEIS A006153 == finite_map_stats.num_depth_le_1.

  4. Asymptotic trend: E[T] grows like sqrt(n) (Flajolet-Odlyzko 1990); we check
     monotonicity of E[T] and of E[T]/sqrt(n) staying in the documented band.

References:
  - Flajolet & Odlyzko (1990), Random Mapping Statistics, EUROCRYPT '89.
  - OEIS A216242 (Geoffrey Critzer), A006153, A000312 (n^n), A000142 (n!).
"""

import math
import sys
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from itertools import product

from src.finite_map_stats import num_depth_le_1
from src.max_tail_length import (
    _max_tail_fast,
    brute_force_distribution,
    closed_form_diagonal_top,
    cumulative_triangle,
    distribution_triangle,
    distribution_via_egf,
    expected_max_tail,
    height_eq_t_counts,
    height_le_t_counts,
    max_tail_length,
    num_permutations,
    tree_egf_height_le,
)


# OEIS A216242: number of functions f:[n]->[n] of height k, read by rows, n>=1.
# (Geoffrey Critzer, Mar 14 2013) -- published DATA, n = 1..9.
A216242_ROWS = [
    [1],
    [2, 2],
    [6, 15, 6],
    [24, 124, 84, 24],
    [120, 1185, 1160, 540, 120],
    [720, 13086, 17610, 10560, 3960, 720],
    [5040, 165361, 296772, 214410, 104160, 32760, 5040],
    [40320, 2363320, 5536440, 4692576, 2686320, 1115520, 302400, 40320],
    [362880, 37780497, 113680800, 111488328, 72080064, 35637840,
     12942720, 3084480, 362880],
]

# OEIS A006153 (height <= 1), offset 1: 1, 4, 21, 148, 1305, 13806, 170401, ...
A006153 = [1, 4, 21, 148, 1305, 13806, 170401, 2403640]

FACTORIAL = [math.factorial(n) for n in range(12)]


# ---------------------------------------------------------------------------
# Single-map primitive
# ---------------------------------------------------------------------------

class TestSingleMap:
    def test_identity_permutation_height_zero(self):
        assert max_tail_length((0, 1, 2, 3)) == 0

    def test_constant_map_height(self):
        # x -> 0 on [4]: 0 is the unique fixed point (cyclic), every other
        # vertex maps directly into it => height 1.
        assert max_tail_length((0, 0, 0, 0)) == 1

    def test_path_into_fixed_point(self):
        # 0->0, 1->0, 2->1, 3->2 : a path 3->2->1->0 with 0 a fixed point.
        # tail lengths: 0:0, 1:1, 2:2, 3:3 => height 3 = n-1.
        assert max_tail_length((0, 0, 1, 2)) == 3

    def test_two_cycle_with_tail(self):
        # 0<->1 a 2-cycle; 2->0, 3->2 : tails 2:1, 3:2 => height 2.
        assert max_tail_length((1, 0, 0, 2)) == 2

    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4, 5, 6])
    def test_fast_path_equals_reference_exhaustively(self, n):
        # brute_force_distribution uses _max_tail_fast; verify it equals the
        # reference max_tail_length on EVERY map for n <= 6.
        for f in product(range(n), repeat=n):
            assert _max_tail_fast(f, n) == max_tail_length(f), (n, f)


# ---------------------------------------------------------------------------
# Core check: EGF route == brute force, term by term
# ---------------------------------------------------------------------------

class TestEgfMatchesBruteForce:
    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7])
    def test_distribution_matches_brute_force(self, n):
        bf = brute_force_distribution(n)
        egf = distribution_via_egf(n)
        assert bf == egf, f"n={n}: brute {bf} != egf {egf}"

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7])
    def test_brute_force_row_equals_triangle_row(self, n):
        bf = brute_force_distribution(n)
        row = distribution_triangle(n)[-1]  # D(n, 0..n-1)
        assert [bf.get(t, 0) for t in range(n)] == row

    @pytest.mark.slow
    def test_distribution_matches_brute_force_n8(self):
        assert brute_force_distribution(8) == distribution_via_egf(8)


# ---------------------------------------------------------------------------
# Boundary / sanity identities
# ---------------------------------------------------------------------------

class TestSanityIdentities:
    @pytest.mark.parametrize("n", range(1, 11))
    def test_T0_equals_n_factorial(self, n):
        # D(n,0) = number of permutations = n!
        assert height_eq_t_counts(0, n)[n] == math.factorial(n)
        assert distribution_via_egf(n)[0] == math.factorial(n)
        assert num_permutations(n) == math.factorial(n)

    @pytest.mark.parametrize("n", range(1, 11))
    def test_top_diagonal_equals_n_factorial(self, n):
        # D(n, n-1) = n!  (a single path of n vertices)
        assert height_eq_t_counts(n - 1, n)[n] == math.factorial(n)
        assert closed_form_diagonal_top(n) == math.factorial(n)

    @pytest.mark.parametrize("n", range(1, 11))
    def test_row_sums_to_n_to_the_n(self, n):
        row = distribution_triangle(n)[-1]
        assert sum(row) == n ** n

    @pytest.mark.parametrize("n", range(1, 11))
    def test_cumulative_saturates_at_n_to_the_n(self, n):
        # H(n, t) = n^n for every t >= n-1.
        for t in range(n - 1, n + 3):
            assert height_le_t_counts(t, n)[n] == n ** n

    @pytest.mark.parametrize("n", range(1, 11))
    def test_cumulative_is_partial_sum_of_exact(self, n):
        row = distribution_triangle(n)[-1]
        cum = cumulative_triangle(n)[-1]
        acc = 0
        for t in range(n):
            acc += row[t]
            assert cum[t] == acc

    def test_height_le_negative_one(self):
        # only the empty map has height <= -1; H(0,-1)=1, H(n,-1)=0 for n>=1.
        col = height_le_t_counts(-1, 6)
        assert col[0] == 1
        assert all(col[n] == 0 for n in range(1, 7))


# ---------------------------------------------------------------------------
# Tree EGF recursion
# ---------------------------------------------------------------------------

class TestTreeEgf:
    def test_e0_is_x(self):
        e0 = tree_egf_height_le(0, 6)
        assert e0[1] == 1
        assert all(e0[k] == 0 for k in range(7) if k != 1)

    def test_tree_counts_height_le_t_are_labeled_rooted_trees(self):
        # As t -> infinity (t >= n-1) the count n! [x^n] e_t must equal the total
        # number of labeled rooted trees on n nodes = n^{n-1} (Cayley).
        N = 8
        e = tree_egf_height_le(N - 1, N)  # height unrestricted up to degree N
        for n in range(1, N + 1):
            c = e[n] * math.factorial(n)
            assert c.denominator == 1
            assert int(c) == n ** (n - 1), f"rooted trees n={n}"

    def test_height_le_1_trees_are_stars(self):
        # rooted trees of height <= 1 on n nodes: root + (n-1) leaves => exactly
        # 1 shape, n labeled (choice of root) => n!*[x^n] e_1 = n.
        N = 7
        e1 = tree_egf_height_le(1, N)
        for n in range(1, N + 1):
            assert int(e1[n] * math.factorial(n)) == n


# ---------------------------------------------------------------------------
# Prior-art / OEIS validation
# ---------------------------------------------------------------------------

class TestOeis:
    def test_distribution_triangle_is_A216242(self):
        rows = distribution_triangle(len(A216242_ROWS))
        assert rows == A216242_ROWS

    def test_height_le_1_is_A006153(self):
        col = [height_le_t_counts(1, n)[n] for n in range(1, len(A006153) + 1)]
        assert col == A006153

    def test_height_le_1_matches_finite_map_stats_depth_le_1(self):
        # num_depth_le_1 in finite_map_stats is the same quantity (A006153).
        for n in range(1, 9):
            assert height_le_t_counts(1, n)[n] == num_depth_le_1(n)

    def test_row_sums_are_n_to_the_n_A000312(self):
        for n, row in enumerate(distribution_triangle(9), start=1):
            assert sum(row) == n ** n


# ---------------------------------------------------------------------------
# Expected max tail length
# ---------------------------------------------------------------------------

class TestExpectation:
    def test_small_exact_values(self):
        # Hand/derived exact rationals (cross-checked against brute force below).
        assert expected_max_tail(1) == Fraction(0)
        assert expected_max_tail(2) == Fraction(1, 2)
        assert expected_max_tail(3) == Fraction(1)        # (0*6+1*15+2*6)/27=27/27
        assert expected_max_tail(4) == Fraction(91, 64)

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7])
    def test_expectation_matches_brute_force(self, n):
        bf = brute_force_distribution(n)
        nn = n ** n
        exp_bf = sum(Fraction(t * c, nn) for t, c in bf.items())
        assert expected_max_tail(n) == exp_bf

    def test_expectation_is_increasing(self):
        vals = [expected_max_tail(n) for n in range(1, 11)]
        assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))

    def test_expectation_grows_like_sqrt_n(self):
        # Flajolet-Odlyzko 1990: E[height] ~ c*sqrt(n).  Check E[T]/sqrt(n) is
        # bounded and increasing toward the documented constant band over a
        # decade of n (slow 1/sqrt(n)-type corrections, so we only bound it).
        ratios = [float(expected_max_tail(n)) / math.sqrt(n) for n in range(4, 11)]
        assert all(0.5 < r < 2.6 for r in ratios)
        assert ratios[-1] > ratios[0]  # still climbing toward the limit


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
