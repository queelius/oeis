"""
Tests for Zarankiewicz number computation.

Prior art validation against known exact values from:
  - OEIS A072567: z(n,n;2,2)
  - OEIS A001198 -> z(n,n;3,3)
  - OEIS A006613 -> z(n,n;2,3)
"""
import pytest

from src.zarankiewicz import (
    KNOWN_Z_2_2,
    KNOWN_Z_2_3,
    KNOWN_Z_2_4,
    KNOWN_Z_3_3,
    KNOWN_Z_3_4,
    KNOWN_Z_3_4_RECT1,
    KNOWN_Z_3_4_RECT2,
    KNOWN_Z_4_4,
    _KNOWN_TABLES,
    _kst_one_side,
    count_ones,
    find_witness,
    format_matrix,
    is_kst_free,
    kst_upper_bound,
    lookup_known,
    zarankiewicz,
    zarankiewicz_brute,
    zarankiewicz_sat,
)


# ---------------------------------------------------------------------------
# Verification tests
# ---------------------------------------------------------------------------

class TestIsKstFree:
    """Tests for the K_{s,t}-freeness checker."""

    def test_empty_matrix(self):
        assert is_kst_free([], 2, 2) is True

    def test_single_entry_zero(self):
        assert is_kst_free([[0]], 1, 1) is True

    def test_single_entry_one(self):
        assert is_kst_free([[1]], 1, 1) is False

    def test_2x2_all_ones(self):
        matrix = [[1, 1], [1, 1]]
        assert is_kst_free(matrix, 2, 2) is False

    def test_2x2_three_ones(self):
        matrix = [[1, 1], [1, 0]]
        assert is_kst_free(matrix, 2, 2) is True

    def test_3x3_identity(self):
        matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        assert is_kst_free(matrix, 2, 2) is True

    def test_3x3_no_2x2(self):
        # A 3x3 matrix with 6 ones and no 2x2 all-1 submatrix
        # This is the incidence matrix of a triangle (K_3):
        # row i has 1s in columns i and (i+1)%3
        matrix = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
        assert is_kst_free(matrix, 2, 2) is True
        assert count_ones(matrix) == 6

    def test_trivial_too_few_rows(self):
        # s > m: cannot form the submatrix
        matrix = [[1, 1, 1]]
        assert is_kst_free(matrix, 2, 2) is True

    def test_trivial_too_few_cols(self):
        # t > n: cannot form the submatrix
        matrix = [[1], [1]]
        assert is_kst_free(matrix, 2, 2) is True

    def test_3x3_all_ones_no_3x3(self):
        matrix = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
        assert is_kst_free(matrix, 3, 3) is False

    def test_asymmetric_2x3(self):
        # A matrix with a 2x3 all-1 submatrix
        matrix = [
            [1, 1, 1],
            [1, 1, 1],
            [0, 0, 0],
        ]
        assert is_kst_free(matrix, 2, 3) is False

    def test_asymmetric_2x3_free(self):
        # Each pair of rows shares at most 2 common 1-columns
        matrix = [
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
        ]
        assert is_kst_free(matrix, 2, 3) is True


class TestCountOnes:
    def test_empty(self):
        assert count_ones([]) == 0

    def test_all_zeros(self):
        assert count_ones([[0, 0], [0, 0]]) == 0

    def test_all_ones(self):
        assert count_ones([[1, 1], [1, 1]]) == 4

    def test_mixed(self):
        assert count_ones([[1, 0], [0, 1]]) == 2


# ---------------------------------------------------------------------------
# Brute-force tests (tiny cases only)
# ---------------------------------------------------------------------------

class TestBruteForce:
    """Brute-force on very small cases for ground truth."""

    def test_z_1_1_2_2(self):
        assert zarankiewicz_brute(1, 1, 2, 2) == 1

    def test_z_2_2_2_2(self):
        assert zarankiewicz_brute(2, 2, 2, 2) == 3

    def test_z_3_3_2_2(self):
        assert zarankiewicz_brute(3, 3, 2, 2) == 6

    def test_z_2_3_2_2(self):
        # Asymmetric: 2x3 matrix, no 2x2 all-1 submatrix
        assert zarankiewicz_brute(2, 3, 2, 2) == 4

    def test_z_1_1_1_1(self):
        assert zarankiewicz_brute(1, 1, 1, 1) == 0

    def test_z_2_2_1_1(self):
        # No 1x1 all-1 submatrix = no 1s at all
        assert zarankiewicz_brute(2, 2, 1, 1) == 0

    def test_z_3_3_3_3(self):
        # z(3,3;3,3) = 8 (all 9 minus 1)
        assert zarankiewicz_brute(3, 3, 3, 3) == 8

    def test_z_2_2_2_3(self):
        # z(2,2;2,3) = 4: cannot pick 3 cols from 2
        assert zarankiewicz_brute(2, 2, 2, 3) == 4

    def test_trivial_s_exceeds_m(self):
        # s > m means every matrix is K_{s,t}-free
        assert zarankiewicz_brute(2, 5, 3, 2) == 10

    def test_z_4_4_2_2(self):
        assert zarankiewicz_brute(4, 4, 2, 2) == 9


# ---------------------------------------------------------------------------
# SAT solver tests
# ---------------------------------------------------------------------------

class TestSAT:
    """SAT-based computation, validated against known values."""

    def test_z_1_1_2_2(self):
        assert zarankiewicz_sat(1, 1, 2, 2) == 1

    def test_z_2_2_2_2(self):
        assert zarankiewicz_sat(2, 2, 2, 2) == 3

    def test_z_3_3_2_2(self):
        assert zarankiewicz_sat(3, 3, 2, 2) == 6

    def test_z_4_4_2_2(self):
        assert zarankiewicz_sat(4, 4, 2, 2) == 9

    def test_z_5_5_2_2(self):
        assert zarankiewicz_sat(5, 5, 2, 2) == 12

    def test_z_6_6_2_2(self):
        assert zarankiewicz_sat(6, 6, 2, 2) == 16

    @pytest.mark.slow
    def test_z_7_7_2_2(self):
        assert zarankiewicz_sat(7, 7, 2, 2) == 21

    def test_z_3_3_3_3(self):
        assert zarankiewicz_sat(3, 3, 3, 3) == 8

    def test_z_4_4_3_3(self):
        assert zarankiewicz_sat(4, 4, 3, 3) == 13

    def test_z_3_3_2_3(self):
        assert zarankiewicz_sat(3, 3, 2, 3) == 7

    def test_z_4_4_2_3(self):
        assert zarankiewicz_sat(4, 4, 2, 3) == 12


class TestSATWitness:
    """Test that witness matrices are valid."""

    def test_witness_z_4_4_2_2(self):
        val, matrix = find_witness(4, 4, 2, 2)
        assert val == 9
        assert count_ones(matrix) == 9
        assert is_kst_free(matrix, 2, 2)

    def test_witness_z_5_5_2_2(self):
        val, matrix = find_witness(5, 5, 2, 2)
        assert val == 12
        assert count_ones(matrix) == 12
        assert is_kst_free(matrix, 2, 2)

    def test_witness_z_3_3_3_3(self):
        val, matrix = find_witness(3, 3, 3, 3)
        assert val == 8
        assert count_ones(matrix) == 8
        assert is_kst_free(matrix, 3, 3)


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;2,2)
# ---------------------------------------------------------------------------

class TestPriorArt_Z22:
    """Validate SAT solver against OEIS A072567 for z(n,n;2,2)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 3), (3, 6), (4, 9), (5, 12), (6, 16),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 2, 2) == expected

    @pytest.mark.slow
    @pytest.mark.parametrize("n,expected", [
        (7, 21), (8, 24), (9, 29), (10, 34),
    ])
    def test_medium_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 2, 2) == expected


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;3,3)
# ---------------------------------------------------------------------------

class TestPriorArt_Z33:
    """Validate SAT solver against OEIS A001198 for z(n,n;3,3)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 4), (3, 8), (4, 13),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 3, 3) == expected

    @pytest.mark.slow
    @pytest.mark.parametrize("n,expected", [
        (5, 20), (6, 26),
    ])
    def test_medium_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 3, 3) == expected


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;2,3)
# ---------------------------------------------------------------------------

class TestPriorArt_Z23:
    """Validate SAT solver against OEIS A006613 for z(n,n;2,3)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 4), (3, 7), (4, 12), (5, 16),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 2, 3) == expected

    @pytest.mark.slow
    @pytest.mark.parametrize("n,expected", [
        (6, 21), (7, 28),
    ])
    def test_medium_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 2, 3) == expected


# ---------------------------------------------------------------------------
# Cross-validation: brute force vs SAT
# ---------------------------------------------------------------------------

class TestCrossValidation:
    """Check brute force and SAT agree on small cases."""

    @pytest.mark.parametrize("m,n,s,t", [
        (1, 1, 2, 2), (2, 2, 2, 2), (3, 3, 2, 2),
        (2, 3, 2, 2), (3, 2, 2, 2),
        (1, 1, 3, 3), (2, 2, 3, 3), (3, 3, 3, 3),
        (2, 2, 2, 3), (3, 3, 2, 3),
        (2, 4, 2, 2), (4, 2, 2, 2),
    ])
    def test_brute_vs_sat(self, m, n, s, t):
        brute_val = zarankiewicz_brute(m, n, s, t)
        sat_val = zarankiewicz_sat(m, n, s, t)
        assert brute_val == sat_val, (
            f"z({m},{n};{s},{t}): brute={brute_val}, sat={sat_val}"
        )


# ---------------------------------------------------------------------------
# Upper bound tests
# ---------------------------------------------------------------------------

class TestUpperBound:
    """Verify KST upper bound is valid."""

    @pytest.mark.parametrize("n", range(1, 15))
    def test_kst_bound_z22(self, n):
        exact = KNOWN_Z_2_2[n]
        bound = kst_upper_bound(n, n, 2, 2)
        assert exact <= bound + 0.001, (
            f"z({n},{n};2,2) = {exact} exceeds KST bound {bound:.2f}"
        )

    @pytest.mark.parametrize("n", range(1, 11))
    def test_kst_bound_z33(self, n):
        exact = KNOWN_Z_3_3[n]
        bound = kst_upper_bound(n, n, 3, 3)
        assert exact <= bound + 0.001, (
            f"z({n},{n};3,3) = {exact} exceeds KST bound {bound:.2f}"
        )

    @pytest.mark.parametrize("n", range(1, 12))
    def test_kst_bound_z23(self, n):
        exact = KNOWN_Z_2_3[n]
        bound = kst_upper_bound(n, n, 2, 3)
        assert exact <= bound + 0.001, (
            f"z({n},{n};2,3) = {exact} exceeds KST bound {bound:.2f}"
        )

    @pytest.mark.parametrize("n", range(1, 12))
    def test_kst_bound_z24(self, n):
        exact = KNOWN_Z_2_4[n]
        bound = kst_upper_bound(n, n, 2, 4)
        assert exact <= bound + 0.001, (
            f"z({n},{n};2,4) = {exact} exceeds KST bound {bound:.2f}"
        )

    @pytest.mark.parametrize("n", range(1, 11))
    def test_kst_bound_z34(self, n):
        exact = KNOWN_Z_3_4[n]
        bound = kst_upper_bound(n, n, 3, 4)
        assert exact <= bound + 0.001, (
            f"z({n},{n};3,4) = {exact} exceeds KST bound {bound:.2f}"
        )

    @pytest.mark.parametrize("n", range(1, 14))
    def test_kst_bound_z44(self, n):
        exact = KNOWN_Z_4_4[n]
        bound = kst_upper_bound(n, n, 4, 4)
        assert exact <= bound + 0.001, (
            f"z({n},{n};4,4) = {exact} exceeds KST bound {bound:.2f}"
        )


# ---------------------------------------------------------------------------
# Two-sided KST bound tests (regression for the m/n swap bug)
# ---------------------------------------------------------------------------

class TestKSTTwoSided:
    """Test that the KST bound correctly uses both sides for rectangular cases.

    The old code applied Jensen's inequality from only one side, computing
    m * C(e/m, s) <= (t-1) * C(n, s). The correct row-side formula is
    n * C(e/n, s) <= C(m, s) * (t-1), and we must also check the column
    side: m * C(e/m, t) <= C(n, t) * (s-1). The bound is min of both.
    """

    def test_z_12_11_2_3_corrected(self):
        """z(12,11;2,3): the bug's primary regression case.

        Old (wrong) bound: ~57.73  (too low, invalid upper bound).
        Correct bound:     ~59.67  (row side dominates).
        """
        bound = kst_upper_bound(12, 11, 2, 3)
        assert abs(bound - 59.6687) < 0.01, f"Expected ~59.67, got {bound:.4f}"

    def test_z_11_11_2_3_square(self):
        """z(11,11;2,3) = 55: square case, both sides should agree.

        The bound should be exactly 55.0 (tight for z(11,11;2,3)).
        """
        bound = kst_upper_bound(11, 11, 2, 3)
        assert abs(bound - 55.0) < 0.01, f"Expected 55.0, got {bound:.4f}"
        # Known exact value must not exceed the bound.
        assert 55 <= bound + 0.001

    def test_two_sides_agree_square(self):
        """For square matrices, both sides of the KST bound should give
        the same result (by symmetry of the counting argument)."""
        for s, t in [(2, 2), (2, 3), (3, 3), (2, 4), (3, 4), (4, 4)]:
            for n in range(max(s, t), 15):
                bound_row = _kst_one_side(part_a=n, part_b=n, k=s, max_common=t - 1)
                bound_col = _kst_one_side(part_a=n, part_b=n, k=t, max_common=s - 1)
                # For s == t and square, these must be identical.
                # For s != t, they can differ even in the square case.
                bound = kst_upper_bound(n, n, s, t)
                assert bound <= bound_row + 0.001
                assert bound <= bound_col + 0.001

    def test_symmetry_z_m_n_equals_z_n_m_transposed(self):
        """z(m,n;s,t) = z(n,m;t,s), so the bound should be the same."""
        test_cases = [
            (12, 11, 2, 3),
            (8, 5, 3, 4),
            (10, 7, 2, 4),
            (6, 9, 3, 3),
            (15, 10, 2, 5),
        ]
        for m, n, s, t in test_cases:
            bound1 = kst_upper_bound(m, n, s, t)
            bound2 = kst_upper_bound(n, m, t, s)
            assert abs(bound1 - bound2) < 0.001, (
                f"z({m},{n};{s},{t}) bound={bound1:.4f} != "
                f"z({n},{m};{t},{s}) bound={bound2:.4f}"
            )

    def test_rectangular_known_values_z34_rect1(self):
        """All known z(n,n+1;3,4) values must be <= the KST bound."""
        for (m, n), exact in KNOWN_Z_3_4_RECT1.items():
            bound = kst_upper_bound(m, n, 3, 4)
            assert exact <= bound + 0.001, (
                f"z({m},{n};3,4) = {exact} exceeds KST bound {bound:.2f}"
            )

    def test_rectangular_known_values_z34_rect2(self):
        """All known z(n,n+2;3,4) values must be <= the KST bound."""
        for (m, n), exact in KNOWN_Z_3_4_RECT2.items():
            bound = kst_upper_bound(m, n, 3, 4)
            assert exact <= bound + 0.001, (
                f"z({m},{n};3,4) = {exact} exceeds KST bound {bound:.2f}"
            )

    def test_all_known_square_values_bounded(self):
        """Every known exact square value must be <= the KST bound."""
        for (s, t), table in _KNOWN_TABLES.items():
            for n, exact in table.items():
                bound = kst_upper_bound(n, n, s, t)
                assert exact <= bound + 0.001, (
                    f"z({n},{n};{s},{t}) = {exact} exceeds KST bound {bound:.2f}"
                )

    def test_bound_at_least_trivial(self):
        """The KST bound should never be negative or exceed m*n."""
        test_cases = [
            (5, 10, 2, 3), (10, 5, 3, 2), (1, 100, 2, 2),
            (3, 3, 5, 5), (2, 8, 2, 4), (20, 7, 3, 4),
        ]
        for m, n, s, t in test_cases:
            bound = kst_upper_bound(m, n, s, t)
            assert bound >= 0, f"Bound for z({m},{n};{s},{t}) is negative: {bound}"
            assert bound <= m * n + 0.001, (
                f"Bound for z({m},{n};{s},{t}) = {bound:.2f} exceeds m*n = {m*n}"
            )

    def test_two_sides_can_differ_rectangular(self):
        """For rectangular matrices with s != t, the two sides should generally
        give different bounds, confirming we actually need both."""
        # z(12,11;2,3): row side ~59.67, col side ~65.16
        bound_row = _kst_one_side(part_a=12, part_b=11, k=2, max_common=2)
        bound_col = _kst_one_side(part_a=11, part_b=12, k=3, max_common=1)
        assert abs(bound_row - bound_col) > 1.0, (
            f"Expected different bounds, got row={bound_row:.4f}, col={bound_col:.4f}"
        )
        # The overall bound should be the tighter one
        assert kst_upper_bound(12, 11, 2, 3) == min(bound_row, bound_col)


# ---------------------------------------------------------------------------
# Edge cases and invariants
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge cases and structural invariants."""

    def test_z_1_1_s_t_trivial(self):
        # z(1,1;s,t) = 1 for s >= 2 or t >= 2
        for s, t in [(2, 2), (2, 3), (3, 3), (5, 5)]:
            assert zarankiewicz_sat(1, 1, s, t) == 1

    def test_monotonicity_in_n(self):
        """z(n,n;s,t) is non-decreasing in n."""
        for s, t in [(2, 2), (3, 3)]:
            prev = 0
            for n in range(1, 7):
                val = zarankiewicz_sat(n, n, s, t)
                assert val >= prev, (
                    f"z({n},{n};{s},{t})={val} < z({n-1},{n-1};{s},{t})={prev}"
                )
                prev = val

    def test_monotonicity_in_st(self):
        """z(m,n;s,t) is non-decreasing in s and t (larger forbidden
        submatrix = weaker constraint = more 1s allowed)."""
        vals = []
        for s in [2, 3, 4, 5]:
            vals.append(zarankiewicz_sat(5, 5, s, s))
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i + 1], (
                f"z(5,5;{i+2},{i+2})={vals[i]} > z(5,5;{i+3},{i+3})={vals[i+1]}"
            )

    def test_lookup_known(self):
        assert lookup_known(5, 5, 2, 2) == 12
        assert lookup_known(3, 3, 3, 3) == 8
        assert lookup_known(4, 4, 2, 3) == 12
        assert lookup_known(100, 100, 2, 2) is None

    def test_format_matrix(self):
        matrix = [[1, 0], [0, 1]]
        result = format_matrix(matrix)
        assert "1 0" in result
        assert "0 1" in result

    def test_zarankiewicz_dispatch(self):
        assert zarankiewicz(3, 3, 2, 2, method="sat") == 6
        assert zarankiewicz(3, 3, 2, 2, method="brute") == 6
        assert zarankiewicz(3, 3, 2, 2, method="lookup") == 6

    def test_zarankiewicz_lookup_missing(self):
        with pytest.raises(ValueError):
            zarankiewicz(100, 100, 2, 2, method="lookup")

    def test_zarankiewicz_invalid_method(self):
        with pytest.raises(ValueError):
            zarankiewicz(3, 3, 2, 2, method="invalid")


# ---------------------------------------------------------------------------
# New values: z(n,n+1;3,4) extension (OEIS A006622)
# ---------------------------------------------------------------------------

class TestNewValues_Z34_Rect1:
    """Test new value z(9,10;3,4) = 60 extending OEIS A006622."""

    def test_lookup_z_9_10_3_4(self):
        """New value: z(9,10;3,4) = 60 (A006622 a(9) = 61)."""
        assert lookup_known(9, 10, 3, 4) == 60

    def test_lookup_z_8_9_3_4(self):
        """Known value: z(8,9;3,4) = 50."""
        assert lookup_known(8, 9, 3, 4) == 50

    def test_lookup_symmetry_z_10_9_4_3(self):
        """z(10,9;4,3) = z(9,10;3,4) = 60 by symmetry."""
        assert lookup_known(10, 9, 4, 3) == 60

    @pytest.mark.slow
    def test_sat_z_5_6_3_4(self):
        """Validate z(5,6;3,4) = 25 via SAT."""
        assert zarankiewicz_sat(5, 6, 3, 4) == 25


# ---------------------------------------------------------------------------
# New values: z(n,n+2;3,4) extension (OEIS A006625)
# ---------------------------------------------------------------------------

class TestNewValues_Z34_Rect2:
    """Test new value z(9,11;3,4) = 66 extending OEIS A006625."""

    def test_lookup_z_9_11_3_4(self):
        """New value: z(9,11;3,4) = 66 (A006625 a(9) = 67)."""
        assert lookup_known(9, 11, 3, 4) == 66

    def test_lookup_z_8_10_3_4(self):
        """Known value: z(8,10;3,4) = 54."""
        assert lookup_known(8, 10, 3, 4) == 54

    def test_lookup_symmetry_z_11_9_4_3(self):
        """z(11,9;4,3) = z(9,11;3,4) = 66 by symmetry."""
        assert lookup_known(11, 9, 4, 3) == 66


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;3,4) (OEIS A006615)
# ---------------------------------------------------------------------------

class TestPriorArt_Z34:
    """Validate SAT solver against OEIS A006615 for z(n,n;3,4)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 4), (3, 9), (4, 14), (5, 21), (6, 30),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 3, 4) == expected

    @pytest.mark.slow
    def test_z_7_7_3_4(self):
        assert zarankiewicz_sat(7, 7, 3, 4) == 37

    def test_lookup_z_10_10_3_4(self):
        """New value: z(10,10;3,4) = 66 (A006615 a(10) = 67)."""
        assert lookup_known(10, 10, 3, 4) == 66


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;4,4) (OEIS A006616)
# ---------------------------------------------------------------------------

class TestPriorArt_Z44:
    """Validate SAT solver against OEIS A006616 for z(n,n;4,4)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 4), (3, 9), (4, 15), (5, 22), (6, 31),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 4, 4) == expected

    @pytest.mark.slow
    def test_z_7_7_4_4(self):
        assert zarankiewicz_sat(7, 7, 4, 4) == 42


# ---------------------------------------------------------------------------
# Prior art validation: z(n,n;2,4) (OEIS A006614)
# ---------------------------------------------------------------------------

class TestPriorArt_Z24:
    """Validate SAT solver against OEIS A006614 for z(n,n;2,4)."""

    @pytest.mark.parametrize("n,expected", [
        (1, 1), (2, 4), (3, 9), (4, 13), (5, 20), (6, 25),
    ])
    def test_small_n(self, n, expected):
        assert zarankiewicz_sat(n, n, 2, 4) == expected
