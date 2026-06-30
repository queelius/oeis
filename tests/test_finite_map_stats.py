"""
Tests for exact statistics of functional graphs of finite maps.

Prior-art validation (>= 3 known quantities, all exact):
  1. total maps           = n^n
  2. number of permutations = n!                       (1,2,6,24,120,720,5040)
  3. idempotent maps      = A000248                    (1,3,10,41,196,1057,6322)
  4. connected functions  = A001865                    (1,3,17,142,1569,21576,355081)
  5. E[# cyclic points]   = sum_k n!/((n-k)! n^k)       (exact Fraction)
  6. sum of fixed points  = n^n                        (mean fixed points = 1)

Restricted-class closed forms (both attested in OEIS, Critzer 2012):
  * depth <= 1 maps       = A006153                    (1,4,21,148,1305,13806,170401)
                            e.g.f. 1/(1 - x e^x); = sum_k C(n,k) k! k^(n-k)
  * in-degree in {0,2}    = A126934 (|a|)              (n=2:2, 4:36, 6:1800, 8:176400)
                            = C(2m,m) (2m)!/2^m

Every closed form is independently checked against brute-force enumeration
over all n^n maps for small n.

References:
  - Flajolet & Odlyzko (1990), Random Mapping Statistics, EUROCRYPT '89.
  - Harris (1960), Probability distributions related to random mappings.
  - OEIS A000248, A001865, A006153, A126934.
"""

import math
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.finite_map_stats import (
    all_maps,
    brute_force_restricted_count,
    cyclic_points,
    enumerate_totals,
    expected_cyclic_points,
    is_depth_le_1,
    is_indegree_in_02,
    map_statistics,
    num_components,
    num_connected,
    num_depth_le_1,
    num_idempotent,
    num_idempotent_by_image_size,
    num_indegree_in_02,
    num_permutations,
    summary_table,
    tail_lengths,
    total_cyclic_points,
    total_maps,
)


# Known OEIS values, offset matched to n = 1, 2, 3, ...
A000248 = [1, 3, 10, 41, 196, 1057, 6322]        # idempotent maps
A001865 = [1, 3, 17, 142, 1569, 21576, 355081]   # connected functions
A006153 = [1, 4, 21, 148, 1305, 13806, 170401]   # depth <= 1 maps (offset 1)
FACTORIAL = [1, 2, 6, 24, 120, 720, 5040]        # permutations = n!


# ---------------------------------------------------------------------------
# Single-map primitives
# ---------------------------------------------------------------------------

class TestSingleMapPrimitives:
    def test_cyclic_points_identity_permutation(self):
        # identity on [4]: every point is a fixed point => all cyclic.
        f = (0, 1, 2, 3)
        assert cyclic_points(f) == {0, 1, 2, 3}

    def test_cyclic_points_constant_map(self):
        # constant map x -> 0 on [4]: only 0 is cyclic (a self-loop fixed pt).
        f = (0, 0, 0, 0)
        assert cyclic_points(f) == {0}

    def test_cyclic_points_pure_cycle(self):
        # 3-cycle 0->1->2->0 on [3]: all cyclic.
        f = (1, 2, 0)
        assert cyclic_points(f) == {0, 1, 2}

    def test_tail_lengths_chain_into_loop(self):
        # 3 -> 2 -> 1 -> 0 -> 0 : 0 is fixed (cyclic), tails 0,1,2,3.
        f = (0, 0, 1, 2)
        assert tail_lengths(f) == [0, 1, 2, 3]

    def test_tail_lengths_all_cyclic(self):
        f = (1, 0)  # transposition, both cyclic
        assert tail_lengths(f) == [0, 0]

    def test_num_components_two_loops(self):
        # 0->0, 1->1 : two components.
        assert num_components((0, 1)) == 2

    def test_num_components_one_tree_into_cycle(self):
        # 0->1->0 cycle, 2->0, 3->2 : single component.
        assert num_components((1, 0, 0, 2)) == 1

    def test_map_statistics_constant(self):
        f = (0, 0, 0, 0)
        st = map_statistics(f)
        assert st["num_cyclic"] == 1
        assert st["num_fixed"] == 1
        assert st["image_size"] == 1
        assert st["num_terminal"] == 3
        assert st["num_components"] == 1
        assert st["max_tail"] == 1
        assert st["is_permutation"] == 0
        assert st["is_idempotent"] == 1     # f(f(x)) = 0 = f(x)
        assert st["is_connected"] == 1

    def test_map_statistics_identity(self):
        f = (0, 1, 2)
        st = map_statistics(f)
        assert st["num_cyclic"] == 3
        assert st["num_fixed"] == 3
        assert st["max_tail"] == 0
        assert st["num_components"] == 3
        assert st["is_permutation"] == 1
        assert st["is_idempotent"] == 1


# ---------------------------------------------------------------------------
# Known quantity 1: total maps = n^n
# ---------------------------------------------------------------------------

class TestTotalMaps:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_total_equals_n_pow_n(self, n):
        assert total_maps(n) == n ** n
        assert enumerate_totals(n)["total_maps"] == n ** n


# ---------------------------------------------------------------------------
# Known quantity 2: permutations = n!
# ---------------------------------------------------------------------------

class TestPermutations:
    @pytest.mark.parametrize("n,expected", list(zip(range(1, 8), FACTORIAL)))
    def test_closed_form_factorial(self, n, expected):
        assert num_permutations(n) == expected

    @pytest.mark.parametrize("n", range(1, 7))
    def test_brute_force_matches_factorial(self, n):
        assert enumerate_totals(n)["counts"]["is_permutation"] == math.factorial(n)


# ---------------------------------------------------------------------------
# Known quantity 3: idempotent maps = A000248
# ---------------------------------------------------------------------------

class TestIdempotent:
    @pytest.mark.parametrize("n,expected", list(zip(range(1, 8), A000248)))
    def test_closed_form_a000248(self, n, expected):
        assert num_idempotent(n) == expected

    @pytest.mark.parametrize("n", range(1, 7))
    def test_brute_force_matches_a000248(self, n):
        assert enumerate_totals(n)["counts"]["is_idempotent"] == num_idempotent(n)

    @pytest.mark.parametrize("n", range(1, 7))
    def test_image_size_refinement_sums(self, n):
        # rows C(n,k) k^(n-k) sum to A000248(n)
        total = sum(num_idempotent_by_image_size(n, k) for k in range(n + 1))
        assert total == num_idempotent(n)

    def test_image_size_refinement_brute_force(self):
        # n=4: idempotent maps by image size must match the closed-form row.
        n = 4
        from collections import Counter
        by_size = Counter()
        for f in product(range(n), repeat=n):
            if all(f[f[x]] == f[x] for x in range(n)):
                by_size[len(set(f))] += 1
        for k in range(n + 1):
            assert by_size[k] == num_idempotent_by_image_size(n, k)


# ---------------------------------------------------------------------------
# Known quantity 4: connected functions = A001865
# ---------------------------------------------------------------------------

class TestConnected:
    @pytest.mark.parametrize("n,expected", list(zip(range(1, 8), A001865)))
    def test_closed_form_a001865(self, n, expected):
        assert num_connected(n) == expected

    @pytest.mark.parametrize("n", range(1, 7))
    def test_brute_force_matches_a001865(self, n):
        assert enumerate_totals(n)["counts"]["is_connected"] == num_connected(n)


# ---------------------------------------------------------------------------
# Known quantity 5: E[# cyclic points] = sum_k n!/((n-k)! n^k)
# ---------------------------------------------------------------------------

class TestExpectedCyclicPoints:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_expectation_is_exact_fraction(self, n):
        e = expected_cyclic_points(n)
        assert isinstance(e, Fraction)
        # match the explicit Q-function sum
        ref = sum(
            Fraction(math.factorial(n), math.factorial(n - k)) * Fraction(1, n ** k)
            for k in range(1, n + 1)
        )
        assert e == ref

    @pytest.mark.parametrize("n", range(1, 7))
    def test_expectation_matches_brute_force(self, n):
        tot = enumerate_totals(n)
        assert tot["expectations"]["num_cyclic"] == expected_cyclic_points(n)

    @pytest.mark.parametrize("n,expected", list(zip(
        range(1, 8), [1, 6, 51, 568, 7845, 129456, 2485567])))
    def test_total_cyclic_points_value(self, n, expected):
        # sum over all maps of (# cyclic points) = n^n * E[# cyclic].
        assert total_cyclic_points(n) == expected
        assert total_cyclic_points(n) == n ** n * expected_cyclic_points(n)

    @pytest.mark.parametrize("n", range(1, 7))
    def test_total_cyclic_points_matches_brute_force(self, n):
        tot = enumerate_totals(n)
        assert total_cyclic_points(n) == tot["sums"]["num_cyclic"]

    @pytest.mark.parametrize("n", range(1, 9))
    def test_total_cyclic_points_is_n_times_connected(self, n):
        # OEIS A063169: a(n) = n * A001865(n).
        assert total_cyclic_points(n) == n * num_connected(n)

    def test_known_small_values(self):
        assert expected_cyclic_points(2) == Fraction(3, 2)
        assert expected_cyclic_points(3) == Fraction(17, 9)
        assert expected_cyclic_points(4) == Fraction(71, 32)


# ---------------------------------------------------------------------------
# Known quantity 6: sum of fixed points = n^n (mean fixed points = 1)
# ---------------------------------------------------------------------------

class TestFixedPoints:
    @pytest.mark.parametrize("n", range(1, 7))
    def test_sum_fixed_equals_n_pow_n(self, n):
        tot = enumerate_totals(n)
        assert tot["sums"]["num_fixed"] == n ** n
        assert tot["expectations"]["num_fixed"] == Fraction(1)


# ---------------------------------------------------------------------------
# Restricted class A: depth <= 1 maps = A006153
# ---------------------------------------------------------------------------

class TestDepthLe1:
    @pytest.mark.parametrize("n,expected", list(zip(range(1, 8), A006153)))
    def test_closed_form_a006153(self, n, expected):
        assert num_depth_le_1(n) == expected

    @pytest.mark.parametrize("n", range(1, 7))
    def test_closed_form_matches_brute_force(self, n):
        assert num_depth_le_1(n) == brute_force_restricted_count(n, is_depth_le_1)

    def test_offset_zero_value(self):
        # A006153 offset 0 has a(0) = 1 (empty map).
        assert num_depth_le_1(0) == 1

    def test_egf_coefficients(self):
        # e.g.f. 1/(1 - x e^x): a(n) = n! [x^n].  Check against sympy series.
        sympy = pytest.importorskip("sympy")
        x = sympy.symbols("x")
        ser = sympy.series(1 / (1 - x * sympy.exp(x)), x, 0, 8).removeO()
        for n in range(1, 8):
            coeff = sympy.factorial(n) * ser.coeff(x, n)
            assert int(coeff) == num_depth_le_1(n)

    def test_depth_le_1_predicate_excludes_depth_2(self):
        # chain 3 -> 2 -> 1 -> 0 -> 0 has a depth-3 vertex; excluded.
        assert not is_depth_le_1((0, 0, 1, 2))
        # but 2 -> 0, 1 -> 0, 0 -> 0 (star into a loop) is all depth <= 1.
        assert is_depth_le_1((0, 0, 0))


# ---------------------------------------------------------------------------
# Restricted class C: in-degree in {0,2} maps = A126934
# ---------------------------------------------------------------------------

class TestIndegreeIn02:
    def test_closed_form_even_n(self):
        # |A126934|: n=2 -> 2, n=4 -> 36, n=6 -> 1800, n=8 -> 176400.
        assert num_indegree_in_02(2) == 2
        assert num_indegree_in_02(4) == 36
        assert num_indegree_in_02(6) == 1800
        assert num_indegree_in_02(8) == 176400

    @pytest.mark.parametrize("n", [1, 3, 5, 7])
    def test_odd_n_is_zero(self, n):
        assert num_indegree_in_02(n) == 0

    @pytest.mark.parametrize("n", range(1, 7))
    def test_closed_form_matches_brute_force(self, n):
        assert num_indegree_in_02(n) == brute_force_restricted_count(n, is_indegree_in_02)

    def test_equals_a001147_times_a001813(self):
        # OEIS A126934: |a(m)| = A001147(m) * A001813(m)
        #   A001147(m) = (2m-1)!!  (double factorial)
        #   A001813(m) = (2m)!/m!  (quadruple-factorial-ish)
        for m in range(1, 5):
            n = 2 * m
            a001147 = math.prod(range(1, 2 * m, 2))          # (2m-1)!!
            a001813 = math.factorial(2 * m) // math.factorial(m)  # (2m)!/m!
            assert num_indegree_in_02(n) == a001147 * a001813


# ---------------------------------------------------------------------------
# Cross-checks and invariants
# ---------------------------------------------------------------------------

class TestInvariants:
    @pytest.mark.parametrize("n", range(1, 6))
    def test_cyclic_plus_tail_equals_n_per_map(self, n):
        # every vertex is cyclic or has positive tail length; counts partition.
        for f in all_maps(n):
            cyc = cyclic_points(f)
            tails = tail_lengths(f)
            non_cyclic = sum(1 for x in range(n) if x not in cyc)
            assert len(cyc) + non_cyclic == n
            # tail length is 0 exactly on cyclic points
            for x in range(n):
                assert (tails[x] == 0) == (x in cyc)

    @pytest.mark.parametrize("n", range(1, 6))
    def test_permutations_have_all_cyclic_no_terminal(self, n):
        for f in all_maps(n):
            if len(set(f)) == n:  # permutation
                assert len(cyclic_points(f)) == n
                assert max(tail_lengths(f)) == 0

    @pytest.mark.parametrize("n", range(2, 6))
    def test_components_le_cyclic_points(self, n):
        # number of components (cycles) <= number of cyclic points.
        for f in all_maps(n):
            assert num_components(f) <= len(cyclic_points(f))

    def test_summary_table_shape(self):
        rows = summary_table(5)
        assert len(rows) == 5
        assert rows[0]["n"] == 1
        # spot-check published columns
        assert rows[4]["connected"] == A001865[4]
        assert rows[4]["depth_le_1"] == A006153[4]
        assert rows[4]["permutations"] == FACTORIAL[4]


@pytest.mark.slow
class TestLargerN:
    def test_n7_connected_and_depth(self):
        # n=7 brute force is ~823k maps; keep behind the slow marker.
        tot = enumerate_totals(7)
        assert tot["counts"]["is_connected"] == A001865[6]
        assert tot["sums"]["num_fixed"] == 7 ** 7
        assert num_depth_le_1(7) == A006153[6]
