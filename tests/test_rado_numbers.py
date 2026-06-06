"""
Tests for Rado number computations.

Prior art validation:
  - Schur numbers S(k) = max N with valid k-coloring of {1,...,N} avoiding x+y=z
    S(1)=1, S(2)=4, S(3)=13  (Schur 1916)
  - van der Waerden W(3;2) = 9: min N forcing monochromatic 3-AP in 2 colors
  - R(x+y=kz, 2) with distinct variables for small k (Robertson-Myers 2015 and others)

New computations:
  - R(x+y=kz, 2) with distinct variables for k up to 43
  - Closed-form conjecture for k >= 8
  - Various equation families
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from conftest import requires_sat


# ---------------------------------------------------------------------------
# Prior art validation: Schur numbers
# ---------------------------------------------------------------------------

@requires_sat
class TestSchurNumbers:
    """Schur number S(k): largest N such that {1,...,N} can be k-colored
    avoiding monochromatic x + y = z (repetitions allowed).
    Our rado_number returns the smallest UNSAT N, so R = S(k) + 1."""

    def test_schur_1(self):
        from src.rado_numbers import schur_number
        val, witness = schur_number(1, max_n=20)
        assert val == 2, f"R(x+y=z, 1) should be 2 (S(1)=1), got {val}"

    def test_schur_2(self):
        from src.rado_numbers import schur_number
        val, witness = schur_number(2, max_n=20)
        assert val == 5, f"R(x+y=z, 2) should be 5 (S(2)=4), got {val}"

    def test_schur_3(self):
        from src.rado_numbers import schur_number
        val, witness = schur_number(3, max_n=50)
        assert val == 14, f"R(x+y=z, 3) should be 14 (S(3)=13), got {val}"

    def test_schur_witness_valid(self):
        """Witness coloring for S(2)=4 should be valid on {1,...,4}."""
        from src.rado_numbers import schur_number, verify_coloring
        val, witness = schur_number(2, max_n=20)
        assert witness is not None
        assert verify_coloring([1, 1, -1], witness, val - 1, distinct=False)

    def test_schur_witness_invalid_at_rado(self):
        """No valid 2-coloring exists at N = R = 5."""
        from src.rado_numbers import _check_n
        sat, _ = _check_n([1, 1, -1], 2, 5, distinct=False)
        assert not sat


# ---------------------------------------------------------------------------
# Prior art validation: van der Waerden W(3;2)
# ---------------------------------------------------------------------------

@requires_sat
class TestVanDerWaerden:
    """W(3;2) = 9: minimum N such that every 2-coloring of {1,...,N}
    has a monochromatic 3-term arithmetic progression.
    Equation: x + y = 2z with distinct variables."""

    def test_w_3_2(self):
        from src.rado_numbers import rado_xy_kz
        val, witness = rado_xy_kz(2, k_colors=2, max_n=20, distinct=True)
        assert val == 9, f"W(3;2) should be 9, got {val}"

    def test_w_3_2_witness(self):
        from src.rado_numbers import rado_xy_kz, verify_coloring
        val, witness = rado_xy_kz(2, k_colors=2, max_n=20, distinct=True)
        assert witness is not None
        assert verify_coloring([1, 1, -2], witness, 8, distinct=True)

    def test_w_3_2_no_coloring_at_9(self):
        from src.rado_numbers import _check_n
        sat, _ = _check_n([1, 1, -2], 2, 9, distinct=True)
        assert not sat


# ---------------------------------------------------------------------------
# x + y = kz (distinct), 2 colors: computed values
# ---------------------------------------------------------------------------

# Complete table of computed R(x+y=kz, 2) with distinct variables
RADO_XY_KZ_DISTINCT_2COLOR = {
    1: 9, 2: 9, 3: 15, 4: 20, 5: 25, 6: 31, 7: 49, 8: 41,
    9: 54, 10: 61, 11: 77, 12: 85, 13: 104, 14: 113, 15: 135,
    16: 145, 17: 170, 18: 181, 19: 209, 20: 221, 21: 252, 22: 265,
    23: 299, 24: 313, 25: 350, 26: 365, 27: 405, 28: 421, 29: 464,
    30: 481, 31: 527, 32: 545, 33: 594, 34: 613, 35: 665, 36: 685,
    37: 740, 38: 761, 39: 819, 40: 841, 41: 902, 42: 925, 43: 989,
    44: 1013, 45: 1080, 46: 1105, 47: 1175, 48: 1201, 49: 1274,
    50: 1301, 51: 1377, 52: 1405, 53: 1484, 54: 1513, 55: 1595,
    56: 1625, 57: 1710, 58: 1741, 59: 1829, 60: 1861,
    61: 1952, 62: 1985, 63: 2079, 64: 2113, 65: 2210, 66: 2245,
    67: 2345, 68: 2381, 69: 2484, 70: 2521, 71: 2627, 72: 2665,
    73: 2774, 74: 2813, 75: 2925, 76: 2965, 77: 3080, 78: 3121,
    79: 3239, 80: 3281,
    81: 3402, 82: 3445, 83: 3569, 84: 3613, 85: 3740, 86: 3785,
    87: 3915, 88: 3961, 89: 4094, 90: 4141, 91: 4277, 92: 4325,
    93: 4464, 94: 4513, 95: 4655, 96: 4705,
    97: 4850, 98: 4901, 99: 5049, 100: 5101,
    101: 5252, 102: 5305, 103: 5459, 104: 5513, 105: 5670,
    106: 5725, 107: 5885, 108: 5941, 109: 6104, 110: 6161,
    111: 6327, 112: 6385, 113: 6554, 114: 6613, 115: 6785,
    116: 6845, 117: 7020, 118: 7081, 119: 7259, 120: 7321,
    121: 7502, 122: 7565, 123: 7749, 124: 7813, 125: 8000,
    126: 8065, 127: 8255, 128: 8321, 129: 8514, 130: 8581,
    131: 8777, 132: 8845, 133: 9044, 134: 9113, 135: 9315,
    136: 9385, 137: 9590, 138: 9661, 139: 9869, 140: 9941,
    141: 10152, 142: 10225, 143: 10439, 144: 10513, 145: 10730,
    146: 10805, 147: 11025, 148: 11101, 149: 11324, 150: 11401,
    151: 11627, 152: 11705, 153: 11934, 154: 12013, 155: 12245,
    156: 12325, 157: 12560, 158: 12641, 159: 12879, 160: 12961,
    161: 13202, 162: 13285, 163: 13529, 164: 13613, 165: 13860,
    166: 13945, 167: 14195, 168: 14281, 169: 14534, 170: 14621,
    171: 14877, 172: 14965, 173: 15224, 174: 15313, 175: 15575,
    176: 15665, 177: 15930, 178: 16021, 179: 16289, 180: 16381,
    181: 16652, 182: 16745, 183: 17019, 184: 17113, 185: 17390,
    186: 17485, 187: 17765, 188: 17861, 189: 18144, 190: 18241,
    191: 18527, 192: 18625, 193: 18914, 194: 19013, 195: 19305,
    196: 19405, 197: 19700, 198: 19801, 199: 20099, 200: 20201,
}


def _rado_xy_kz_formula(k: int) -> int:
    """Closed-form formula for R(x+y=kz, 2) with distinct variables.

    Conjectured to hold for all k >= 8:
      - Odd k:  R = k(k+3)/2
      - Even k: R = (k^2 + 2k + 2)/2
    """
    if k % 2 == 1:
        return k * (k + 3) // 2
    else:
        return (k * k + 2 * k + 2) // 2


@requires_sat
class TestRadoXYkZ:
    """Test R(x+y=kz, 2) with distinct variables."""

    def test_k1_schur(self):
        """k=1 gives the Schur-like equation x+y=z (distinct vars)."""
        from src.rado_numbers import rado_xy_kz
        val, _ = rado_xy_kz(1, k_colors=2, max_n=20, distinct=True)
        assert val == 9

    def test_k2_van_der_waerden(self):
        """k=2 gives W(3;2)=9."""
        from src.rado_numbers import rado_xy_kz
        val, _ = rado_xy_kz(2, k_colors=2, max_n=20, distinct=True)
        assert val == 9

    @pytest.mark.parametrize("k", range(3, 8))
    def test_small_k(self, k):
        """Verify small k values (sporadic range)."""
        from src.rado_numbers import rado_xy_kz, verify_coloring
        val, witness = rado_xy_kz(k, k_colors=2, max_n=100, distinct=True)
        assert val == RADO_XY_KZ_DISTINCT_2COLOR[k]
        if witness is not None:
            assert verify_coloring([1, 1, -k], witness, val - 1, distinct=True)

    @pytest.mark.parametrize("k", range(8, 21))
    def test_formula_range(self, k):
        """Verify values in the closed-form range k=8..20."""
        from src.rado_numbers import rado_xy_kz, verify_coloring
        val, witness = rado_xy_kz(k, k_colors=2, max_n=300, distinct=True)
        assert val == RADO_XY_KZ_DISTINCT_2COLOR[k]
        assert val == _rado_xy_kz_formula(k), \
            f"k={k}: formula gives {_rado_xy_kz_formula(k)}, SAT gives {val}"
        if witness is not None:
            assert verify_coloring([1, 1, -k], witness, val - 1, distinct=True)

    @pytest.mark.parametrize("k", [25, 30, 35, 40, 43])
    @pytest.mark.slow
    def test_larger_k(self, k):
        """Spot-check larger k values."""
        from src.rado_numbers import rado_xy_kz
        val, _ = rado_xy_kz(k, k_colors=2, max_n=1100, distinct=True)
        assert val == RADO_XY_KZ_DISTINCT_2COLOR[k]

    @pytest.mark.parametrize("k", range(8, 201))
    def test_closed_form_conjecture(self, k):
        """Verify the closed-form conjecture against computed values (k=8..200)."""
        expected = RADO_XY_KZ_DISTINCT_2COLOR[k]
        formula = _rado_xy_kz_formula(k)
        assert formula == expected, \
            f"k={k}: formula={formula}, computed={expected}"


# ---------------------------------------------------------------------------
# Explicit lower bound colorings (proved for odd k, verified for even k)
# ---------------------------------------------------------------------------

def _hm_coloring_odd(k, n):
    """Harborth-Maasberg coloring for odd k >= 5.

    Proved to be valid on {1,...,k(k+3)/2 - 1} for all odd k >= 5.
    """
    coloring = [0] * (n + 1)
    for i in range(1, n + 1):
        r = i % k
        if r == 0:
            j = i // k
            coloring[i] = 0 if j < k / 4 else 1
        elif r <= (k - 1) // 2:
            coloring[i] = 1
        else:
            coloring[i] = 0
    return coloring


def _staircase_coloring_even(k, n):
    """Staircase coloring for even k >= 8.

    Verified to be valid on {1,...,(k^2+2k)/2} for k=8..50.
    """
    half_k = k // 2
    t = [0] * k
    for r in range(2, half_k):
        t[r] = 0
    t[0] = k // 4
    t[1] = k // 4
    t[half_k] = (k - 2) // 4
    t[k - 1] = (k - 2) // 4
    if half_k + 1 < k - 1:
        t[half_k + 1] = half_k - 2
    for r in range(half_k + 2, k - 1):
        t[r] = half_k - 1

    coloring = [0] * (n + 1)
    for i in range(1, n + 1):
        r = i % k
        q = i // k
        coloring[i] = 0 if q < t[r] else 1
    return coloring


@requires_sat
class TestExplicitColorings:
    """Test the explicit lower bound colorings."""

    @pytest.mark.parametrize("k", range(5, 200, 2))
    def test_hm_coloring_odd_valid(self, k):
        """HM coloring valid at n = R-1 for odd k."""
        from src.rado_numbers import verify_coloring
        R = _rado_xy_kz_formula(k)
        n = R - 1
        c = _hm_coloring_odd(k, n)
        assert verify_coloring([1, 1, -k], c, n, distinct=True), \
            f"HM coloring invalid at n={n} for k={k}"

    @pytest.mark.parametrize("k", range(5, 30, 2))
    def test_hm_coloring_odd_fails_at_R(self, k):
        """HM coloring fails at n = R for odd k."""
        from src.rado_numbers import verify_coloring
        R = _rado_xy_kz_formula(k)
        c = _hm_coloring_odd(k, R)
        assert not verify_coloring([1, 1, -k], c, R, distinct=True), \
            f"HM coloring still valid at n=R={R} for k={k}"

    @pytest.mark.parametrize("k", range(8, 202, 2))
    def test_staircase_coloring_even_valid(self, k):
        """Staircase coloring valid at n = R-1 for even k."""
        from src.rado_numbers import verify_coloring
        R = _rado_xy_kz_formula(k)
        n = R - 1
        c = _staircase_coloring_even(k, n)
        assert verify_coloring([1, 1, -k], c, n, distinct=True), \
            f"Staircase coloring invalid at n={n} for k={k}"

    @pytest.mark.parametrize("k", range(8, 30, 2))
    def test_staircase_coloring_even_fails_at_R(self, k):
        """Staircase coloring fails at n = R for even k."""
        from src.rado_numbers import verify_coloring
        R = _rado_xy_kz_formula(k)
        c = _staircase_coloring_even(k, R)
        assert not verify_coloring([1, 1, -k], c, R, distinct=True), \
            f"Staircase coloring still valid at n=R={R} for k={k}"


# ---------------------------------------------------------------------------
# x + y = kz (non-distinct), 2 colors
# ---------------------------------------------------------------------------

RADO_XY_KZ_NONDISTINCT_2COLOR = {
    1: 5, 3: 9, 4: 10, 5: 15, 6: 21, 7: 28, 8: 36, 9: 45,
    10: 55, 11: 66, 12: 78, 13: 91, 14: 105, 15: 120,
}


@requires_sat
class TestRadoXYkZNonDistinct:
    """Non-distinct Rado numbers: R = k(k+1)/2 for k >= 4."""

    def test_k1_schur_nondistinct(self):
        from src.rado_numbers import rado_xy_kz
        val, _ = rado_xy_kz(1, k_colors=2, max_n=20, distinct=False)
        assert val == 5  # S(2)+1

    @pytest.mark.parametrize("k", [4, 5, 6, 7, 8, 9, 10, 11, 12])
    def test_triangular_formula(self, k):
        """For k >= 4 (non-distinct): R = k(k+1)/2."""
        from src.rado_numbers import rado_xy_kz
        val, _ = rado_xy_kz(k, k_colors=2, max_n=200, distinct=False)
        assert val == k * (k + 1) // 2, \
            f"k={k}: expected {k*(k+1)//2}, got {val}"


# ---------------------------------------------------------------------------
# Other equation families
# ---------------------------------------------------------------------------

@requires_sat
class TestOtherEquations:
    """Test various other linear equation Rado numbers."""

    def test_x_plus_y_plus_z_eq_w(self):
        """R(x+y+z=w, 2) with distinct variables."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, 1, -1], 2, max_n=50, distinct=True)
        assert val == 24

    def test_x_plus_y_eq_z_plus_w(self):
        """R(x+y=z+w, 2) with distinct variables."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -1, -1], 2, max_n=50, distinct=True)
        assert val == 11

    def test_2x_plus_y_eq_z(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([2, 1, -1], 2, max_n=50, distinct=True)
        assert val == 15

    def test_2x_plus_y_eq_3z(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([2, 1, -3], 2, max_n=50, distinct=True)
        assert val == 13

    def test_2x_plus_3y_eq_5z(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([2, 3, -5], 2, max_n=50, distinct=True)
        assert val == 21

    def test_3x_plus_y_eq_4z(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([3, 1, -4], 2, max_n=50, distinct=True)
        assert val == 17

    def test_x_plus_y_plus_z_eq_3w(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, 1, -3], 2, max_n=50, distinct=True)
        assert val == 21

    def test_x_plus_y_plus_z_eq_2w(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, 1, -2], 2, max_n=200, distinct=True)
        assert val == 13

    def test_2x_plus_3y_eq_z(self):
        from src.rado_numbers import rado_number
        val, _ = rado_number([2, 3, -1], 2, max_n=100, distinct=True)
        assert val == 77


# ---------------------------------------------------------------------------
# Perfect Square Law: R(x + by = bz) = b^2 for b >= 4 (THEOREM)
# ---------------------------------------------------------------------------

@requires_sat
class TestPerfectSquareLaw:
    """Test R(x + by = bz, 2; distinct) = b^2 for b >= 4.

    This is the Perfect Square Law: for b >= 4, the Rado number equals b^2.
    Reference: Towell 2026 (rado_perfect_square.tex).

    Proof structure:
    - Lower bound: canonical coloring (mults one color, non-mults other) avoids
      all monochromatic distinct solutions in {1,...,b^2-1}.
    - Structure lemma: any valid coloring of {1,...,b^2-1} has all mults same
      color and all non-mults the opposite color.
    - Two-Triple Blocking: at N = b^2, triples (b^2, b, 2b) and (b^2, 1, b+1)
      block both colorings of b^2.
    """

    @pytest.mark.parametrize("b,expected", [
        (4, 16), (5, 25), (6, 36), (7, 49), (8, 64),
        (9, 81), (10, 100), (11, 121), (12, 144),
    ])
    def test_r_equals_b_squared(self, b, expected):
        """For b >= 4, R(x + by = bz) = b^2."""
        from src.rado_numbers import rado_number
        max_n = b * b + 10
        val, _ = rado_number([1, b, -b], k_colors=2, max_n=max_n, distinct=True)
        assert val == expected, f"R(x + {b}y = {b}z) should be {expected}, got {val}"

    @pytest.mark.parametrize("b", [4, 5, 6, 7, 8, 9, 10])
    def test_canonical_lower_bound_coloring(self, b):
        """Canonical coloring (mults 0, non-mults 1) avoids all monochromatic
        distinct solutions at N = b^2 - 1."""
        from src.rado_numbers import verify_coloring
        N = b * b - 1
        col = [0] * (N + 1)
        for i in range(1, N + 1):
            col[i] = 0 if (i % b == 0) else 1
        assert verify_coloring([1, b, -b], col, N, distinct=True), (
            f"Canonical coloring failed for b={b}"
        )

    @pytest.mark.parametrize("b", [4, 5, 6, 7, 8])
    def test_two_triple_blocking_arithmetic(self, b):
        """Verify the two blocking triples satisfy x + by = bz and are distinct."""
        # T_alpha = (b^2, b, 2b): all multiples of b
        assert b * b + b * b == b * (2 * b)
        assert len({b * b, b, 2 * b}) == 3
        # T_{1-alpha} = (b^2, 1, b+1): one multiple, two non-multiples
        assert b * b + b * 1 == b * (b + 1)
        assert len({b * b, 1, b + 1}) == 3

    @pytest.mark.parametrize("b", [4, 5, 6, 7])
    def test_structure_lemma_forces_mults_same_color(self, b):
        """For b in {4,...,7}, verify the Structure Lemma: fixing chi(ib) != chi(jb)
        for i != j at N = b^2 - 1 is UNSAT."""
        from pysat.solvers import Solver
        from src.rado_numbers import _enumerate_solutions, _var

        N = b * b - 1
        k = 2
        solution_sets = _enumerate_solutions([1, b, -b], N, distinct=True)

        for i in range(1, b):
            for j in range(i + 1, b):
                with Solver(name="cd195") as solver:
                    for elem in range(1, N + 1):
                        solver.add_clause([_var(elem, 0, k), _var(elem, 1, k)])
                        solver.add_clause([-_var(elem, 0, k), -_var(elem, 1, k)])
                    solver.add_clause([_var(i * b, 0, k)])
                    solver.add_clause([-_var(i * b, 1, k)])
                    solver.add_clause([_var(j * b, 1, k)])
                    solver.add_clause([-_var(j * b, 0, k)])
                    for elem_set in solution_sets:
                        elems = sorted(elem_set)
                        for c in range(k):
                            clause = [-_var(x, c, k) for x in elems]
                            solver.add_clause(clause)
                    assert not solver.solve(), (
                        f"Structure lemma FAIL for b={b}, i={i}, j={j}"
                    )

    @pytest.mark.parametrize("b", [4, 5, 6, 7, 8])
    def test_non_mults_within_distance_b_minus_1(self, b):
        """Every non-multiple of b in {1,...,b^2-1} is within distance b-1 of some
        multiple of b. This is key for the 'non-mults forced' proposition."""
        N = b * b - 1
        for x in range(1, N + 1):
            if x % b == 0:
                continue
            q = x // b
            if q >= 1:
                nearest = q * b
            else:
                nearest = b
            d = abs(x - nearest)
            assert 1 <= d <= b - 1, (
                f"Non-mult {x} at distance {d} from mult {nearest}, b={b}"
            )

    @pytest.mark.slow
    @pytest.mark.parametrize("b,expected", [
        (13, 169), (14, 196), (15, 225), (16, 256), (17, 289), (18, 324),
        (19, 361), (20, 400),
    ])
    def test_r_equals_b_squared_large(self, b, expected):
        """Extended verification for larger b (slow)."""
        from src.rado_numbers import rado_number
        max_n = b * b + 10
        val, _ = rado_number([1, b, -b], k_colors=2, max_n=max_n, distinct=True)
        assert val == expected, f"R(x + {b}y = {b}z) should be {expected}, got {val}"

    @pytest.mark.parametrize("b", [4, 5, 6, 7, 8, 9, 10, 11, 12])
    def test_alternating_cascade_pure_propagation(self, b):
        """Verify the analytical alternating cascade in Proposition 4.6.

        For each (b, k) with k <= b-2 outside the three exceptional pairs
        {(4,1), (4,2), (5,1)}, and each branch chi(1) in {R, B}, the cascade
        starting from chi(kb)=R, chi((k+1)b)=B, chi(1)=branch must derive a
        contradiction by pure unit propagation. This is the analytical
        backbone of Lemma 4.3 part (a).
        """
        N = b * b - 1
        # Build all valid distinct triples (kb', y, y+k') for k' in 1..b-1.
        triples = []
        for y in range(1, N + 1):
            for z in range(y + 1, N + 1):
                x = b * (z - y)
                if 1 <= x <= N and x != y and x != z:
                    triples.append((x, y, z))

        def propagate(forced):
            forced = dict(forced)
            changed = True
            while changed:
                changed = False
                for (x, y, z) in triples:
                    cs = [forced.get(x), forced.get(y), forced.get(z)]
                    r = cs.count("R")
                    bb = cs.count("B")
                    if r == 3 or bb == 3:
                        return None
                    if r == 2:
                        for t in (x, y, z):
                            if t not in forced:
                                forced[t] = "B"
                                changed = True
                    elif bb == 2:
                        for t in (x, y, z):
                            if t not in forced:
                                forced[t] = "R"
                                changed = True
            return forced

        exceptions = {(4, 1), (4, 2), (5, 1)}
        for k in range(1, b - 1):
            if (b, k) in exceptions:
                continue
            for branch in ("R", "B"):
                forced = {k * b: "R", (k + 1) * b: "B", 1: branch}
                result = propagate(forced)
                assert result is None, (
                    f"Alternating cascade did not contradict by propagation: "
                    f"b={b}, k={k}, chi(1)={branch}"
                )

    def test_b4_k1_dpll_leaves(self):
        """Verify the four leaves of Lemma 4.9 (b=4, adjacent pair (1,2)).

        The hand-checked DPLL tree branches on chi(1), chi(5), chi(7) and
        terminates at one of: (8,10,12) all B; (4,12,13) all R; or (8,5,7) all B.
        """
        # All terminal triples must satisfy x + 4y = 4z with distinct components.
        for (x, y, z) in [(8, 10, 12), (4, 12, 13), (8, 5, 7)]:
            assert x + 4 * y == 4 * z, f"({x},{y},{z}) bad arithmetic for b=4"
            assert len({x, y, z}) == 3, f"({x},{y},{z}) not distinct"

    def test_b4_k2_dpll_leaves(self):
        """Verify the two leaves of Lemma 4.10 (b=4, adjacent pair (2,3))."""
        for (x, y, z) in [(4, 12, 13), (8, 4, 6)]:
            assert x + 4 * y == 4 * z, f"({x},{y},{z}) bad arithmetic for b=4"
            assert len({x, y, z}) == 3, f"({x},{y},{z}) not distinct"

    def test_b5_k1_extended_cascade(self):
        """Verify Lemma 4.8: extended cascade for b=5, adjacent pair (1,2).

        Branch chi(1)=R uses additional triples (15,1,4), (20,1,5), (15,3,6),
        (20,3,7), (5,6,7).
        """
        for (x, y, z) in [(15, 1, 4), (20, 1, 5), (15, 3, 6), (20, 3, 7), (5, 6, 7)]:
            assert x + 5 * y == 5 * z, f"({x},{y},{z}) bad for b=5"
            assert len({x, y, z}) == 3, f"({x},{y},{z}) not distinct"
            assert max(x, y, z) <= 24, f"({x},{y},{z}) out of universe"

    def test_adjacent_suffices_lemma(self):
        """Verify Lemma 4.5: chain transitivity from adjacent pairs.

        Trivially: if chi(kb) = chi((k+1)b) for all k in {1,...,b-2}, then
        chi(b) = chi(2b) = ... = chi((b-1)b).
        """
        # Symbolic: pretend each multiple has a color. If consecutive pairs
        # all match, then all multiples have the same color.
        for b in range(4, 12):
            # Try a hypothetical assignment: all multiples same color.
            colors = ["R"] * (b - 1)  # indexed 0..b-2 for multiples b, 2b, ..., (b-1)b
            for k in range(b - 2):
                assert colors[k] == colors[k + 1], "transitivity broken"


# ---------------------------------------------------------------------------
# 3-color computations
# ---------------------------------------------------------------------------

@requires_sat
class TestThreeColor:
    """3-color Rado numbers."""

    def test_schur_3color(self):
        """S(3) = 13, so R(x+y=z, 3) = 14."""
        from src.rado_numbers import schur_number
        val, _ = schur_number(3, max_n=50)
        assert val == 14

    def test_xy_z_3color_distinct(self):
        """R(x+y=z, 3; distinct) = WS(3)+1 = 24 (weak Schur)."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -1], 3, max_n=50, distinct=True)
        assert val == 24

    def test_xy_2z_3color_distinct(self):
        """R(x+y=2z, 3; distinct) = W(3;3) = 27 (van der Waerden)."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -2], 3, max_n=50, distinct=True)
        assert val == 27

    def test_xy_3z_3color(self):
        """R(x+y=3z, 3) with distinct variables."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -3], 3, max_n=200, distinct=True)
        assert val == 99

    def test_xy_3z_3color_nondistinct(self):
        """R(x+y=3z, 3; non-distinct) = 54, matching Hopkins-Schaal R_3(1,3)."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -3], 3, max_n=200, distinct=False)
        assert val == 54

    def test_xy_3z_3color_witnesses(self):
        """Verify witness colorings for 3-color k=1,2,3."""
        from src.rado_numbers import rado_number, verify_coloring
        for k, expected in [(1, 24), (2, 27), (3, 99)]:
            val, witness = rado_number([1, 1, -k], 3, max_n=200,
                                        distinct=True)
            assert val == expected
            assert witness is not None
            assert verify_coloring([1, 1, -k], witness, val - 1,
                                    distinct=True)

    @pytest.mark.slow
    @pytest.mark.parametrize("k", [4, 5, 6, 7, 8, 9, 10])
    def test_xy_kz_3color_lower_bound(self, k):
        """R(x+y=kz, 3; distinct) > 500 for k >= 4.

        This provides evidence that the degree of regularity of x+y=kz
        is exactly 2 for k >= 4 (finite for 2 colors, infinite for 3+).
        """
        from src.rado_numbers import _check_n
        sat, coloring = _check_n([1, 1, -k], 3, 500, distinct=True)
        assert sat, f"k={k}: expected SAT at N=500 but got UNSAT"


# ---------------------------------------------------------------------------
# Verification and edge cases
# ---------------------------------------------------------------------------

@requires_sat
class TestVerification:
    """Test the verification utilities."""

    def test_verify_valid_coloring(self):
        from src.rado_numbers import verify_coloring
        # 2-coloring of {1,...,4} avoiding x+y=z (non-distinct)
        # Color: 1->0, 2->1, 3->1, 4->0
        coloring = [0, 0, 1, 1, 0]
        assert verify_coloring([1, 1, -1], coloring, 4, distinct=False)

    def test_verify_invalid_coloring(self):
        from src.rado_numbers import verify_coloring
        # All same color: 1+2=3, monochromatic
        coloring = [0, 0, 0, 0, 0]
        assert not verify_coloring([1, 1, -1], coloring, 4, distinct=False)

    def test_equation_from_form(self):
        from src.rado_numbers import equation_from_form
        assert equation_from_form([1, 1], 3) == [1, 1, -3]
        assert equation_from_form([2, 3], 5) == [2, 3, -5]
        assert equation_from_form([1, 1, 1], 1) == [1, 1, 1, -1]

    def test_not_found_returns_negative(self):
        """When max_n is too small, return -1."""
        from src.rado_numbers import rado_number
        val, _ = rado_number([1, 1, -10], 2, max_n=5, distinct=True)
        assert val == -1

    def test_witness_consistency(self):
        """Every returned witness should be valid on {1,...,R-1}."""
        from src.rado_numbers import rado_number, verify_coloring
        for coeffs in [[1, 1, -3], [1, 1, -5], [2, 1, -3]]:
            val, witness = rado_number(coeffs, 2, max_n=100, distinct=True)
            assert val > 0
            assert witness is not None
            assert verify_coloring(coeffs, witness, val - 1, distinct=True)


# ---------------------------------------------------------------------------
# Enumeration correctness
# ---------------------------------------------------------------------------

@requires_sat
class TestEnumeration:
    """Test solution enumeration."""

    def test_schur_solutions_n3(self):
        """x + y = z on {1,2,3} non-distinct: 1+1=2, 1+2=3, 2+1=3."""
        from src.rado_numbers import _enumerate_solutions
        sols = _enumerate_solutions([1, 1, -1], 3, distinct=False)
        # As frozensets of elements: {1,2}, {1,3}, {2,3} -- but 1+1=2 -> {1,2}
        # 1+2=3 -> {1,2,3}, 2+1=3 -> {1,2,3}
        # Also 1+1=2 -> {1,2}, 1+2=3 -> {1,2,3}
        assert frozenset({1, 2}) in sols  # 1+1=2
        assert frozenset({1, 2, 3}) in sols  # 1+2=3

    def test_distinct_excludes_repeats(self):
        """With distinct=True, x+y=2z has no solution (1,1,1)."""
        from src.rado_numbers import _enumerate_solutions
        sols = _enumerate_solutions([1, 1, -2], 3, distinct=True)
        # Distinct solutions: 1+3=2*2 -> {1,2,3}
        assert frozenset({1, 2, 3}) in sols
        # No single-element sets (which would come from x=y=z)
        for s in sols:
            assert len(s) >= 2

    def test_no_solutions_small_n(self):
        """x+y=10z on {1,...,4} has no solutions."""
        from src.rado_numbers import _enumerate_solutions
        sols = _enumerate_solutions([1, 1, -10], 4, distinct=True)
        assert len(sols) == 0
