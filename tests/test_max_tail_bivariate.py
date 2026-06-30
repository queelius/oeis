"""
Tests for the BIVARIATE refinements of the max-tail-length distribution of
functional graphs:

    B(n,t,c)  = #{f : max tail length t, exactly c cyclic points}
    Bm(n,t,m) = #{f : max tail length t, exactly m components}

Validation strategy (exact arithmetic throughout):

  1. JOINT vs BRUTE FORCE.  Each bivariate EGF is checked cell-by-cell against
     exhaustive enumeration over all n^n maps, n = 1..6 (n = 7 under `slow`).
       (t,c):  EGF  1/(1 - u e_t)
       (t,m):  EGF  exp(u log(1/(1-e_t))) = (1 - e_t)^{-u}

  2. MARGINALS, both validated exactly:
       sum_c B(n,t,c) = D(n,t)            (A216242, univariate max tail length)
       sum_t B(n,t,c) = M(n,c)            (A066324, endofunctions by cyclic pts)
       sum_m Bm(n,t,m)= D(n,t)            (A216242)
       sum_t Bm(n,t,m)= A060281(n,m)      (endofunctions by components / cycles)

  3. STRUCTURAL identities (proved in the bivariate README):
       row sums  sum_{t,c} B(n,t,c) = n^n
       t=0 column is the permutation column: B(n,0,c) = [c=n] n!  (one cyclic
            class: all n points cyclic);  Bm(n,0,m) = Stirling1 |s(n,m)| (a
            permutation with m cycles), summing to n!.
       support: B(n,t,c)=0 unless c>=1 and (t=0 ⟹ c=n) and (t>=1 ⟹ c<=n-t).

  4. PRIOR-ART / OEIS: the two marginal triangles equal their published OEIS
     DATA exactly (A066324, A060281); the univariate marginal equals A216242.

References:
  - Flajolet & Sedgewick, "Analytic Combinatorics" (2009), II.4, III.7.
  - OEIS A216242 (Critzer), A066324 (cyclic points), A060281 (components),
    A000312 (n^n), A008275/A130534 (Stirling numbers of the first kind).
"""

import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.max_tail_bivariate import (
    brute_force_joint_tc,
    brute_force_joint_tm,
    cum_tc_egf,
    cum_tm_egf,
    cyclic_points_marginal,
    flatten_joint_tc_support,
    flatten_joint_tm_support,
    joint_tc_rows,
    joint_tc_via_egf,
    joint_tm_rows,
    joint_tm_via_egf,
)
from src.max_tail_length import distribution_triangle


# ---------------------------------------------------------------------------
# Published OEIS DATA for the two marginal triangles (rows n = 1..)
# ---------------------------------------------------------------------------

# A066324: endofunctions on [n] with exactly c (= k) cyclic points (rooted
# trees in the cycle), read by rows.
A066324_ROWS = [
    [1],
    [2, 2],
    [9, 12, 6],
    [64, 96, 72, 24],
    [625, 1000, 900, 480, 120],
    [7776, 12960, 12960, 8640, 3600, 720],
    [117649, 201684, 216090, 164640, 88200, 30240, 5040],
]

# A060281: endofunctions on [n] with exactly m (= k) components / cycles.
A060281_ROWS = [
    [1],
    [3, 1],
    [17, 9, 1],
    [142, 95, 18, 1],
    [1569, 1220, 305, 30, 1],
    [21576, 18694, 5595, 745, 45, 1],
    [355081, 334369, 113974, 18515, 1540, 63, 1],
]

# A216242: univariate max-tail-length triangle (cross-check the marginal).
A216242_ROWS = [
    [1],
    [2, 2],
    [6, 15, 6],
    [24, 124, 84, 24],
    [120, 1185, 1160, 540, 120],
    [720, 13086, 17610, 10560, 3960, 720],
    [5040, 165361, 296772, 214410, 104160, 32760, 5040],
]

# |Stirling1(n,m)| (unsigned): permutations of [n] with m cycles -- the t=0
# row of the (t,m) joint. A130534 read by rows, n = 1..
STIRLING1_ABS_ROWS = [
    [1],
    [1, 1],
    [2, 3, 1],
    [6, 11, 6, 1],
    [24, 50, 35, 10, 1],
    [120, 274, 225, 85, 15, 1],
    [720, 1764, 1624, 735, 175, 21, 1],
]


# ---------------------------------------------------------------------------
# 1. Joint EGF == brute force, cell by cell
# ---------------------------------------------------------------------------

class TestJointTcMatchesBruteForce:
    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
    def test_tc_egf_equals_brute_force(self, n):
        bf = brute_force_joint_tc(n)
        eg = joint_tc_via_egf(n)
        assert bf == eg, f"n={n}: brute {bf} != egf {eg}"

    @pytest.mark.slow
    def test_tc_egf_equals_brute_force_n7(self):
        assert brute_force_joint_tc(7) == joint_tc_via_egf(7)


class TestJointTmMatchesBruteForce:
    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
    def test_tm_egf_equals_brute_force(self, n):
        bf = brute_force_joint_tm(n)
        eg = joint_tm_via_egf(n)
        assert bf == eg, f"n={n}: brute {bf} != egf {eg}"

    @pytest.mark.slow
    def test_tm_egf_equals_brute_force_n7(self):
        assert brute_force_joint_tm(7) == joint_tm_via_egf(7)


# ---------------------------------------------------------------------------
# 2. Marginals -- both directions, both joints, exact
# ---------------------------------------------------------------------------

class TestMarginals:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_tc_marginal_over_c_is_A216242(self, n):
        # sum_c B(n,t,c) = D(n,t), the univariate max-tail distribution.
        j = joint_tc_via_egf(n)
        marg = Counter()
        for (t, c), v in j.items():
            marg[t] += v
        assert [marg[t] for t in range(n)] == distribution_triangle(n)[-1]
        assert [marg[t] for t in range(n)] == A216242_ROWS[n - 1]

    @pytest.mark.parametrize("n", range(1, 8))
    def test_tc_marginal_over_t_is_A066324(self, n):
        # sum_t B(n,t,c) = M(n,c), endofunctions by number of cyclic points.
        j = joint_tc_via_egf(n)
        marg = Counter()
        for (t, c), v in j.items():
            marg[c] += v
        row = [marg[c] for c in range(1, n + 1)]
        assert row == A066324_ROWS[n - 1]
        assert row == cyclic_points_marginal(n)  # closed form

    @pytest.mark.parametrize("n", range(1, 8))
    def test_tm_marginal_over_m_is_A216242(self, n):
        j = joint_tm_via_egf(n)
        marg = Counter()
        for (t, m), v in j.items():
            marg[t] += v
        assert [marg[t] for t in range(n)] == A216242_ROWS[n - 1]

    @pytest.mark.parametrize("n", range(1, 8))
    def test_tm_marginal_over_t_is_A060281(self, n):
        j = joint_tm_via_egf(n)
        marg = Counter()
        for (t, m), v in j.items():
            marg[m] += v
        assert [marg[m] for m in range(1, n + 1)] == A060281_ROWS[n - 1]


# ---------------------------------------------------------------------------
# 3. Structural identities
# ---------------------------------------------------------------------------

class TestStructuralIdentities:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_total_is_n_to_the_n_tc(self, n):
        assert sum(joint_tc_via_egf(n).values()) == n ** n

    @pytest.mark.parametrize("n", range(1, 8))
    def test_total_is_n_to_the_n_tm(self, n):
        assert sum(joint_tm_via_egf(n).values()) == n ** n

    @pytest.mark.parametrize("n", range(1, 8))
    def test_t0_column_tc_is_permutation_all_cyclic(self, n):
        # T = 0 iff permutation iff every point cyclic: B(n,0,c) = [c=n] n!.
        j = joint_tc_via_egf(n)
        for c in range(1, n + 1):
            expected = math.factorial(n) if c == n else 0
            assert j.get((0, c), 0) == expected

    @pytest.mark.parametrize("n", range(1, 8))
    def test_t0_row_tm_is_unsigned_stirling1(self, n):
        # T = 0 maps are exactly the permutations; by number of cycles this is
        # |Stirling1(n,m)|, summing to n!.
        j = joint_tm_via_egf(n)
        row = [j.get((0, m), 0) for m in range(1, n + 1)]
        assert row == STIRLING1_ABS_ROWS[n - 1]
        assert sum(row) == math.factorial(n)

    @pytest.mark.parametrize("n", range(1, 8))
    def test_tc_support_is_constrained(self, n):
        # Nonzero only when c >= 1, and (t == 0 => c == n), (t >= 1 => c <= n-t).
        for (t, c), v in joint_tc_via_egf(n).items():
            assert v > 0
            assert 1 <= c <= n
            if t == 0:
                assert c == n
            else:
                assert c <= n - t

    @pytest.mark.parametrize("n", range(2, 8))
    def test_top_tail_diagonal_is_single_path(self, n):
        # t = n-1 forces a single path of n vertices ending in a fixed point:
        # one cyclic point, one component, n! labelings.
        jtc = joint_tc_via_egf(n)
        jtm = joint_tm_via_egf(n)
        assert jtc[(n - 1, 1)] == math.factorial(n)
        assert sum(v for (t, c), v in jtc.items() if t == n - 1) == math.factorial(n)
        assert jtm[(n - 1, 1)] == math.factorial(n)


# ---------------------------------------------------------------------------
# 4. Closed-form marginal helper and its OEIS identity
# ---------------------------------------------------------------------------

class TestCyclicPointsMarginalClosedForm:
    @pytest.mark.parametrize("n", range(1, 9))
    def test_closed_form_matches_A066324_and_sums_to_n_to_n(self, n):
        row = cyclic_points_marginal(n)
        assert sum(row) == n ** n
        if n <= len(A066324_ROWS):
            assert row == A066324_ROWS[n - 1]

    def test_top_marginal_is_n_factorial(self):
        # c = n cyclic points => permutation => n! maps.
        for n in range(1, 9):
            assert cyclic_points_marginal(n)[-1] == math.factorial(n)


# ---------------------------------------------------------------------------
# 5. Bivariate EGF internal consistency (cumulative -> diff -> nonneg)
# ---------------------------------------------------------------------------

class TestEgfInternals:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_cum_tc_saturates_to_univariate_cyclic_marginal(self, n):
        # At t = n-1 (height unrestricted) the cumulative (t,c) EGF gives the
        # full cyclic-points marginal M(n,c).
        nfac = math.factorial(n)
        cum = cum_tc_egf(n - 1, n)[n]
        row = [int(cum.get(c, Fraction(0)) * nfac) for c in range(1, n + 1)]
        assert row == cyclic_points_marginal(n)

    @pytest.mark.parametrize("n", range(1, 8))
    def test_cum_tm_saturates_to_components_marginal(self, n):
        nfac = math.factorial(n)
        cum = cum_tm_egf(n - 1, n)[n]
        row = [int(cum.get(m, Fraction(0)) * nfac) for m in range(1, n + 1)]
        assert row == A060281_ROWS[n - 1]

    @pytest.mark.parametrize("n", range(1, 7))
    def test_all_joint_counts_nonnegative(self, n):
        assert all(v >= 0 for v in joint_tc_via_egf(n).values())
        assert all(v >= 0 for v in joint_tm_via_egf(n).values())


# ---------------------------------------------------------------------------
# 6. Triangle / flatten plumbing
# ---------------------------------------------------------------------------

class TestTrianglesAndFlatten:
    def test_joint_tc_rows_support_shape(self):
        # n=4 support-only rows: t=0 -> [n!], then c=1..n-t.
        rows = joint_tc_rows(4)
        assert rows[0] == [24]                 # t=0, c=4 only
        assert rows[1] == [4, 48, 72]          # t=1, c=1..3
        assert rows[2] == [36, 48]             # t=2, c=1..2
        assert rows[3] == [24]                 # t=3, c=1

    def test_joint_tc_rows_full_square_pads_zero(self):
        rows = joint_tc_rows(4, full_square=True)
        assert rows[0] == [0, 0, 0, 24]
        assert rows[1] == [4, 48, 72, 0]
        assert all(len(r) == 4 for r in rows)

    def test_joint_tm_rows_t0_is_stirling1(self):
        rows = joint_tm_rows(5)
        assert rows[0] == STIRLING1_ABS_ROWS[4]   # t=0 row, m=1..5

    def test_joint_tm_rows_full_square_pads_zero(self):
        rows = joint_tm_rows(5, full_square=True)
        assert all(len(r) == 5 for r in rows)
        assert rows[0] == STIRLING1_ABS_ROWS[4]   # t=0 already full width
        assert rows[-1] == [120, 0, 0, 0, 0]      # t=4: single path, m=1

    def test_flatten_tc_support_prefix(self):
        flat = flatten_joint_tc_support(4)
        # n=1: [1]; n=2: [2],[2]; n=3: [6],[3,12],[6]; n=4: [24],[4,48,72],[36,48],[24]
        assert flat == [1, 2, 2, 6, 3, 12, 6, 24, 4, 48, 72, 36, 48, 24]

    def test_flatten_tc_support_total_equals_sum_powers(self):
        # The flattened support contains every nonzero count exactly once, so it
        # sums to sum_{n=1}^{N} n^n.
        N = 6
        flat = flatten_joint_tc_support(N)
        assert sum(flat) == sum(n ** n for n in range(1, N + 1))

    def test_flatten_tm_support_total_equals_sum_powers(self):
        N = 6
        flat = flatten_joint_tm_support(N)
        assert sum(flat) == sum(n ** n for n in range(1, N + 1))

    def test_summary_runs(self, capsys):
        from src.max_tail_bivariate import summary
        summary(5)
        out = capsys.readouterr().out
        assert "A216242" in out and "A066324" in out
        assert "True" in out and "False" not in out  # all marginal checks pass


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
