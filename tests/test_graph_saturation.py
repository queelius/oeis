"""
Tests for graph saturation number computation.

Prior art validation:
  - sat(n, K_t) = (t-2)*n - C(t-1, 2)  (Erdos-Hajnal-Moon 1964)
  - sat(n, C_4) = floor((3n-5)/2)  (Ollmann 1972)
  - sat(n, P_4) = n/2 if n even, (n+3)/2 if n odd  (Kaszonyi-Tuza 1986)
  - sat(n, C_5) = ceil(10*(n-1)/7) with exceptions  (Chen 2009/2011)
"""

import pytest
import networkx as nx

from src.graph_saturation import (
    complete_graph,
    contains_subgraph,
    cycle_graph,
    is_saturated,
    path_graph,
    sat_formula_c4,
    sat_formula_c5,
    sat_formula_complete,
    sat_formula_p4,
    sat_number,
)


# ---------------------------------------------------------------------------
# Helper checks
# ---------------------------------------------------------------------------

class TestContainsSubgraph:
    """Test the subgraph containment checker."""

    def test_triangle_in_k4(self):
        assert contains_subgraph(nx.complete_graph(4), nx.complete_graph(3))

    def test_no_triangle_in_tree(self):
        assert not contains_subgraph(nx.star_graph(5), nx.complete_graph(3))

    def test_c5_in_k5(self):
        assert contains_subgraph(nx.complete_graph(5), nx.cycle_graph(5))

    def test_c5_in_k5_minus_edge(self):
        """K5 minus one edge still contains C5."""
        G = nx.complete_graph(5)
        G.remove_edge(0, 1)
        assert contains_subgraph(G, nx.cycle_graph(5))

    def test_no_c5_in_bipartite(self):
        """Bipartite graphs have no odd cycles."""
        G = nx.complete_bipartite_graph(3, 3)
        assert not contains_subgraph(G, nx.cycle_graph(5))

    def test_c4_in_k4(self):
        assert contains_subgraph(nx.complete_graph(4), nx.cycle_graph(4))

    def test_no_c4_in_tree(self):
        assert not contains_subgraph(nx.star_graph(10), nx.cycle_graph(4))

    def test_p4_in_path(self):
        assert contains_subgraph(nx.path_graph(5), nx.path_graph(4))

    def test_no_p4_in_star(self):
        """Star has no path of length 3."""
        assert not contains_subgraph(nx.star_graph(5), nx.path_graph(4))

    def test_empty_graph(self):
        G = nx.Graph()
        G.add_nodes_from(range(5))
        assert not contains_subgraph(G, nx.complete_graph(2))

    def test_self_containment(self):
        G = nx.cycle_graph(5)
        assert contains_subgraph(G, nx.cycle_graph(5))

    def test_non_induced_c5(self):
        """C5 plus a chord still contains C5 as non-induced subgraph."""
        G = nx.cycle_graph(5)
        G.add_edge(0, 2)  # chord
        assert contains_subgraph(G, nx.cycle_graph(5))


class TestIsSaturated:
    """Test the saturation checker."""

    def test_star_is_k3_saturated(self):
        """Star K_{1,n-1} is K_3-saturated."""
        for n in range(3, 8):
            G = nx.star_graph(n - 1)
            assert is_saturated(G, complete_graph(3)), f"Star on {n} vertices"

    def test_complete_graph_not_saturated(self):
        """K_n is not K_3-saturated (it contains K_3)."""
        assert not is_saturated(nx.complete_graph(5), complete_graph(3))

    def test_empty_graph_not_saturated(self):
        """Empty graph on 5 vertices is not K_3-saturated."""
        G = nx.Graph()
        G.add_nodes_from(range(5))
        assert not is_saturated(G, complete_graph(3))

    def test_too_few_vertices_not_complete(self):
        """A non-complete graph with fewer vertices than F is NOT
        F-saturated because adding a non-edge doesn't create F."""
        G = nx.Graph()
        G.add_nodes_from(range(3))
        G.add_edge(0, 1)
        # Adding edge (0,2) or (1,2) doesn't create K4
        assert not is_saturated(G, complete_graph(4))

    def test_vacuous_saturation(self):
        """K_n with n < |V(F)| is vacuously F-saturated
        (F-free and no non-edges to check)."""
        G = nx.complete_graph(2)
        # K_2 is K_3-free and has no non-edges: vacuously saturated
        assert is_saturated(G, complete_graph(3))


# ---------------------------------------------------------------------------
# Known formulas
# ---------------------------------------------------------------------------

class TestFormulas:
    """Test that formula functions compute known values."""

    def test_complete_formula(self):
        """Erdos-Hajnal-Moon: sat(n, K_t) = (t-2)*n - C(t-1,2)."""
        assert sat_formula_complete(5, 3) == 4  # n-1 for triangles
        assert sat_formula_complete(10, 3) == 9
        assert sat_formula_complete(5, 4) == 7
        assert sat_formula_complete(10, 4) == 17
        # t=5, n=5: (5-2)*5 - C(4,2) = 15-6 = 9
        # K5-{e} is K5-saturated on 5 vertices (9 edges).
        assert sat_formula_complete(5, 5) == 9

    def test_c4_formula(self):
        """Ollmann: sat(n, C_4) = floor((3n-5)/2)."""
        assert sat_formula_c4(5) == 5
        assert sat_formula_c4(6) == 6  # (18-5)//2 = 6
        assert sat_formula_c4(7) == 8
        assert sat_formula_c4(10) == 12  # (30-5)//2 = 12
        assert sat_formula_c4(4) is None

    def test_p4_formula(self):
        """Kaszonyi-Tuza: sat(n, P_4)."""
        assert sat_formula_p4(4) == 2
        assert sat_formula_p4(5) == 4  # (5+3)/2 = 4
        assert sat_formula_p4(6) == 3  # 6/2 = 3
        assert sat_formula_p4(7) == 5  # (7+3)/2 = 5
        assert sat_formula_p4(8) == 4  # 8/2 = 4
        assert sat_formula_p4(3) is None

    def test_c5_formula(self):
        """Chen: sat(n, C_5) = ceil(10*(n-1)/7) with exceptions."""
        assert sat_formula_c5(5) == 6
        assert sat_formula_c5(6) == 8
        # n=11 is exceptional: ceil(100/7) - 1 = 15 - 1 = 14
        assert sat_formula_c5(11) == 14
        # n=21 is not exceptional: ceil(200/7) = 29
        assert sat_formula_c5(21) == 29
        assert sat_formula_c5(4) is None


# ---------------------------------------------------------------------------
# SAT solver validation against known formulas
# ---------------------------------------------------------------------------

class TestSatNumberComplete:
    """Verify sat(n, K_t) against Erdos-Hajnal-Moon formula."""

    @pytest.mark.parametrize("n", range(3, 8))
    def test_k3(self, n):
        """sat(n, K_3) = n - 1."""
        result = sat_number(n, complete_graph(3))
        expected = sat_formula_complete(n, 3)
        assert result == expected, f"sat({n}, K3) = {result}, expected {expected}"

    @pytest.mark.parametrize("n", range(4, 9))
    def test_k4(self, n):
        """sat(n, K_4) = 2n - 3."""
        result = sat_number(n, complete_graph(4))
        expected = sat_formula_complete(n, 4)
        assert result == expected, f"sat({n}, K4) = {result}, expected {expected}"

    def test_k4_returns_saturated_graph(self):
        """The returned graph should actually be K_4-saturated."""
        n = 7
        result, G = sat_number(n, complete_graph(4), return_graph=True)
        assert result == sat_formula_complete(n, 4)
        assert is_saturated(G, complete_graph(4))


class TestSatNumberC4:
    """Verify sat(n, C_4) against Ollmann formula."""

    @pytest.mark.parametrize("n", range(5, 9))
    def test_c4_formula(self, n):
        """sat(n, C_4) = floor((3n-5)/2)."""
        result = sat_number(n, cycle_graph(4))
        expected = sat_formula_c4(n)
        assert result == expected, f"sat({n}, C4) = {result}, expected {expected}"

    @pytest.mark.slow
    @pytest.mark.parametrize("n", [9, 10])
    def test_c4_formula_large(self, n):
        """sat(n, C_4) for larger n (slower)."""
        result = sat_number(n, cycle_graph(4))
        expected = sat_formula_c4(n)
        assert result == expected, f"sat({n}, C4) = {result}, expected {expected}"


class TestSatNumberC5:
    """Verify sat(n, C_5) against Chen formula."""

    @pytest.mark.parametrize("n", range(5, 8))
    def test_c5_formula(self, n):
        """sat(n, C_5) matches Chen's formula for small n."""
        result = sat_number(n, cycle_graph(5))
        expected = sat_formula_c5(n)
        assert result == expected, f"sat({n}, C5) = {result}, expected {expected}"

    def test_c5_returns_saturated_graph(self):
        """The returned graph should actually be C_5-saturated."""
        n = 6
        result, G = sat_number(n, cycle_graph(5), return_graph=True)
        assert is_saturated(G, cycle_graph(5))

    @pytest.mark.slow
    @pytest.mark.parametrize("n", [8, 9])
    def test_c5_formula_large(self, n):
        """sat(n, C_5) for larger n (slower)."""
        result = sat_number(n, cycle_graph(5))
        expected = sat_formula_c5(n)
        assert result == expected, f"sat({n}, C5) = {result}, expected {expected}"


class TestSatNumberP4:
    """Verify sat(n, P_4) against Kaszonyi-Tuza formula."""

    @pytest.mark.parametrize("n", range(4, 10))
    def test_p4_formula(self, n):
        """sat(n, P_4) matches Kaszonyi-Tuza formula."""
        result = sat_number(n, path_graph(4))
        expected = sat_formula_p4(n)
        assert result == expected, f"sat({n}, P4) = {result}, expected {expected}"


class TestSatNumberC6:
    """Compute sat(n, C_6) -- no known exact formula.

    Asymptotic bounds (Lan-Shi 2021): 4n/3 - 2 <= sat(n, C6) <= (4n+1)/3 for n >= 9.
    """

    def test_c6_n6(self):
        """sat(6, C6) = 9."""
        result, G = sat_number(6, cycle_graph(6), return_graph=True)
        assert result == 9, f"sat(6, C6) = {result}, expected 9"
        assert is_saturated(G, cycle_graph(6))

    @pytest.mark.slow
    def test_c6_n7(self):
        """sat(7, C6) = 10 (slower computation)."""
        result, G = sat_number(7, cycle_graph(6), return_graph=True)
        assert result == 10, f"sat(7, C6) = {result}, expected 10"
        assert is_saturated(G, cycle_graph(6))


class TestSatNumberP5:
    """Compute sat(n, P_5) for small n."""

    @pytest.mark.parametrize("n,expected", [(5, 4), (6, 5), (7, 6)])
    def test_p5_computed_values(self, n, expected):
        """Computationally verified values for small n."""
        result, G = sat_number(n, path_graph(5), return_graph=True)
        assert result == expected, f"sat({n}, P5) = {result}, expected {expected}"
        assert is_saturated(G, path_graph(5))

    @pytest.mark.slow
    def test_p5_n8(self):
        """sat(8, P5) = 6 (same as n=7; slower computation)."""
        result, G = sat_number(8, path_graph(5), return_graph=True)
        assert result == 6
        assert is_saturated(G, path_graph(5))


class TestSatNumberK23:
    """Compute sat(n, K_{2,3}) for small n."""

    @pytest.mark.parametrize("n,expected", [(5, 7), (6, 9)])
    def test_k23_computed_values(self, n, expected):
        """Computationally verified values for small n."""
        K23 = nx.complete_bipartite_graph(2, 3)
        result, G = sat_number(n, K23, return_graph=True)
        assert result == expected, f"sat({n}, K_{{2,3}}) = {result}, expected {expected}"
        assert is_saturated(G, K23)

    @pytest.mark.slow
    def test_k23_n7(self):
        """sat(7, K_{2,3}) = 11 (slower computation)."""
        K23 = nx.complete_bipartite_graph(2, 3)
        result, G = sat_number(7, K23, return_graph=True)
        assert result == 11
        assert is_saturated(G, K23)


# ---------------------------------------------------------------------------
# Edge cases and invariants
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge cases and structural invariants."""

    def test_n_less_than_f(self):
        """sat(n, F) = None when n < |V(F)|."""
        assert sat_number(2, complete_graph(3)) is None
        assert sat_number(3, complete_graph(4)) is None
        assert sat_number(3, cycle_graph(4)) is None
        assert sat_number(4, cycle_graph(5)) is None

    def test_n_equals_f(self):
        """sat(n, F) when n = |V(F)| -- the minimum case."""
        # sat(3, K3): K3-saturated on 3 vertices.
        # The only K3-free graph on 3 vertices with adding any edge creating K3
        # is the graph with 2 edges (path P3).
        result = sat_number(3, complete_graph(3))
        assert result == sat_formula_complete(3, 3)  # 1*3 - C(2,2) = 2
        assert result == 2

    def test_returned_graph_verified(self):
        """Every returned graph must pass is_saturated."""
        test_cases = [
            (5, complete_graph(3)),
            (6, complete_graph(4)),
            (6, cycle_graph(4)),
            (5, path_graph(4)),
            (6, cycle_graph(5)),
        ]
        for n, F in test_cases:
            result, G = sat_number(n, F, return_graph=True)
            assert is_saturated(G, F), \
                f"Graph returned for n={n} not F-saturated"

    def test_saturation_number_monotonicity(self):
        """sat(n+1, F) >= sat(n, F) is NOT always true for sat numbers,
        but sat(n, K_t) is monotonically increasing."""
        prev = None
        for n in range(3, 8):
            val = sat_number(n, complete_graph(3))
            if prev is not None:
                assert val >= prev, \
                    f"sat({n}, K3) = {val} < sat({n-1}, K3) = {prev}"
            prev = val
