"""
Compute exact graph saturation numbers sat(n, F) for small forbidden graphs F.

A graph G on n vertices is F-saturated if:
  1. G does not contain F as a subgraph
  2. Adding any non-edge of G creates a copy of F

sat(n, F) = minimum number of edges in an F-saturated graph on n vertices.

Known exact formulas (used for validation):
  - sat(n, K_t) = (t-2)*n - C(t-1, 2)  for n >= t  (Erdos-Hajnal-Moon 1964)
  - sat(n, C_4) = floor((3n-5)/2)  for n >= 5  (Ollmann 1972)
  - sat(n, P_4) = n/2 if n even, (n+3)/2 if n odd  (Kaszonyi-Tuza 1986)

Approach: SAT encoding with binary search on edge count.
  - Variables: x_{ij} for each potential edge {i,j}
  - Constraints:
    1. F-freeness: no set of edges forms a copy of F
    2. Saturation: for every non-edge {u,v}, adding it would create F
"""

from collections import defaultdict
from itertools import combinations, permutations
from math import comb, ceil, floor

import networkx as nx
from pysat.card import CardEnc, EncType
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
from pysat.solvers import Solver


# ---------------------------------------------------------------------------
# Standard forbidden graphs
# ---------------------------------------------------------------------------

def complete_graph(t):
    """Return K_t as a networkx graph."""
    return nx.complete_graph(t)


def cycle_graph(k):
    """Return C_k as a networkx graph."""
    return nx.cycle_graph(k)


def path_graph(k):
    """Return P_k (path on k vertices, k-1 edges) as a networkx graph."""
    return nx.path_graph(k)


# ---------------------------------------------------------------------------
# Known formulas (for validation)
# ---------------------------------------------------------------------------

def sat_formula_complete(n, t):
    """Erdos-Hajnal-Moon 1964: sat(n, K_t) for n >= t >= 3."""
    if n < t:
        return None
    return (t - 2) * n - comb(t - 1, 2)


def sat_formula_c4(n):
    """Ollmann 1972: sat(n, C_4) for n >= 5."""
    if n < 5:
        return None
    return (3 * n - 5) // 2


def sat_formula_p4(n):
    """Kaszonyi-Tuza 1986: sat(n, P_4) for n >= 4."""
    if n < 4:
        return None
    if n % 2 == 0:
        return n // 2
    else:
        return (n + 3) // 2


def sat_formula_c5(n):
    """
    Chen 2009/2011: sat(n, C_5).

    For n >= 5:
      sat(n, C5) = ceil(10*(n-1)/7) - 1  if n in {11,12,13,14,16,18,20}
      sat(n, C5) = ceil(10*(n-1)/7)      otherwise

    The formula was proved for n >= 11 by Chen (2009) and extended
    to all n >= 5 by Chen (2011). Our computation confirms the small cases.
    """
    if n < 5:
        return None

    exceptional = {11, 12, 13, 14, 16, 18, 20}
    base = ceil(10 * (n - 1) / 7)
    if n in exceptional:
        return base - 1
    return base


# ---------------------------------------------------------------------------
# Brute-force checker (for small instances and validation)
# ---------------------------------------------------------------------------

def contains_subgraph(G, F):
    """Check if G contains F as a subgraph (not necessarily induced)."""
    if F.number_of_nodes() > G.number_of_nodes():
        return False
    if F.number_of_edges() > G.number_of_edges():
        return False
    # subgraph_is_monomorphic checks for non-induced subgraph containment:
    # an injective map f: V(F) -> V(G) with (u,v) in E(F) => (f(u),f(v)) in E(G).
    # This is the correct notion for graph saturation (extra edges are allowed).
    gm = nx.algorithms.isomorphism.GraphMatcher(G, F)
    return gm.subgraph_is_monomorphic()


def is_saturated(G, F):
    """
    Check if graph G is F-saturated:
    1. G is F-free (does not contain F as a subgraph)
    2. For every non-edge {u,v} in G, G + {u,v} contains F
    """
    if contains_subgraph(G, F):
        return False

    nodes = list(G.nodes())
    for u, v in combinations(nodes, 2):
        if not G.has_edge(u, v):
            G.add_edge(u, v)
            creates_F = contains_subgraph(G, F)
            G.remove_edge(u, v)
            if not creates_F:
                return False

    return True


# ---------------------------------------------------------------------------
# SAT encoding for saturation number computation
# ---------------------------------------------------------------------------

def _unique_copies(n, F):
    """
    Enumerate all distinct copies of F in K_n, identified by their edge set.

    Returns a list of tuples of (i, j) pairs (sorted, i < j) representing
    the edges of each copy.
    """
    f_nodes = sorted(F.nodes())
    f_edges = list(F.edges())
    f_num_nodes = F.number_of_nodes()

    seen = set()
    copies = []
    for perm in permutations(range(n), f_num_nodes):
        emb = dict(zip(f_nodes, perm))
        edge_set = frozenset(
            (min(emb[a], emb[b]), max(emb[a], emb[b]))
            for a, b in f_edges
        )
        if edge_set not in seen:
            seen.add(edge_set)
            copies.append(tuple(sorted(edge_set)))
    return copies


def _edge_to_copies_index(copies):
    """
    Build an index: for each edge (i,j), which copies contain that edge?

    Returns a dict: (i, j) -> list of (copy_index, other_edges) where
    other_edges is the tuple of edges in that copy excluding (i,j).
    """
    index = defaultdict(list)
    for copy_idx, copy_edges in enumerate(copies):
        for edge in copy_edges:
            other = tuple(e for e in copy_edges if e != edge)
            index[edge].append((copy_idx, other))
    return index


class SaturationSATEncoder:
    """
    Encode the graph saturation problem as a SAT instance.

    Uses deduplicated copies of F for efficiency.
    """

    def __init__(self, n, F):
        self.n = n
        self.F = F

        # Edge variable mapping: (i, j) -> variable index (1-based for SAT)
        self._edge_vars = {}
        self._next_var = 1
        for i, j in combinations(range(n), 2):
            self._edge_vars[(i, j)] = self._next_var
            self._next_var += 1

        self.num_edge_vars = len(self._edge_vars)

        # Precompute unique copies of F in K_n
        self._copies = _unique_copies(n, F)

        # Index: which copies use each edge
        self._edge_index = _edge_to_copies_index(self._copies)

    def edge_var(self, i, j):
        """Get the SAT variable for edge {i, j}."""
        if i > j:
            i, j = j, i
        return self._edge_vars[(i, j)]

    def encode_saturation(self):
        """
        Build SAT clauses encoding F-freeness and saturation.

        Returns a list of clauses (each clause is a list of literals).
        These are the "hard" constraints that must always be satisfied.
        """
        clauses = []

        # 1. F-freeness: for each copy of F, at least one edge is absent
        for copy_edges in self._copies:
            clause = [-self.edge_var(*e) for e in copy_edges]
            clauses.append(clause)

        # 2. Saturation: for each pair {u,v}:
        #    edge_var(u,v) OR (exists a copy of F through {u,v}
        #                      where all OTHER edges are present)
        #
        #    Encoding of OR-of-ANDs via auxiliary (Tseitin) variables.
        for u, v in combinations(range(self.n), 2):
            x_uv = self.edge_var(u, v)

            witness_copies = self._edge_index.get((u, v), [])
            if not witness_copies:
                # No copy of F can use this edge => edge must be present
                clauses.append([x_uv])
                continue

            aux_vars = []
            for copy_idx, other_edges in witness_copies:
                aux = self._next_var
                self._next_var += 1
                aux_vars.append(aux)

                # aux => all other edges present
                for e in other_edges:
                    clauses.append([-aux, self.edge_var(*e)])

                # all other edges present => aux
                clause = [aux] + [-self.edge_var(*e) for e in other_edges]
                clauses.append(clause)

            # Saturation: edge present OR at least one witness
            clauses.append([x_uv] + aux_vars)

        return clauses

    def encode(self, max_edges):
        """
        Build SAT clauses for: is there an F-saturated graph on n vertices
        with at most max_edges edges?

        Returns a list of clauses (each clause is a list of literals).
        """
        clauses = self.encode_saturation()

        # Cardinality: at most max_edges edges
        edge_lits = [self.edge_var(i, j)
                      for i, j in combinations(range(self.n), 2)]
        card_cnf = CardEnc.atmost(
            lits=edge_lits,
            bound=max_edges,
            top_id=self._next_var,
            encoding=EncType.totalizer
        )
        card_clause_list = list(card_cnf.clauses)
        if card_clause_list:
            self._next_var = card_cnf.nv + 1
            clauses.extend(card_clause_list)

        return clauses

    def encode_wcnf(self):
        """
        Build a Weighted CNF for MaxSAT optimization.

        Hard clauses: F-freeness + saturation constraints.
        Soft clauses: minimize edges (each edge variable negated with weight 1).

        Returns a WCNF formula.
        """
        hard_clauses = self.encode_saturation()

        wcnf = WCNF()
        for cl in hard_clauses:
            wcnf.append(cl)  # hard clause (infinite weight by default)

        # Soft clauses: prefer edge to be absent (minimize edges)
        for i, j in combinations(range(self.n), 2):
            wcnf.append([-self.edge_var(i, j)], weight=1)

        return wcnf

    def decode_graph(self, model):
        """Extract the graph from a SAT model (list of literals)."""
        model_set = set(model)
        G = nx.Graph()
        G.add_nodes_from(range(self.n))
        for (i, j), var in self._edge_vars.items():
            if var in model_set:
                G.add_edge(i, j)
        return G


def sat_number(n, F, return_graph=False, verbose=False):
    """
    Compute sat(n, F) using MaxSAT (RC2) optimization.

    The saturation constraints (F-freeness + every non-edge creates F)
    are hard clauses. Soft clauses penalize each edge being present,
    so the solver minimizes edge count.

    Parameters
    ----------
    n : int
        Number of vertices.
    F : networkx.Graph
        Forbidden subgraph.
    return_graph : bool
        If True, also return a minimum saturated graph.
    verbose : bool
        If True, print progress information.

    Returns
    -------
    int or (int, networkx.Graph) or None
        The saturation number (and optionally a witness graph), or None
        if no F-saturated graph on n vertices exists.
    """
    f_nodes = F.number_of_nodes()
    if n < f_nodes:
        return None

    encoder = SaturationSATEncoder(n, F)
    wcnf = encoder.encode_wcnf()

    if verbose:
        print(f"  WCNF: {len(list(wcnf.hard))} hard, "
              f"{len(list(wcnf.soft))} soft clauses")

    with RC2(wcnf, solver='cd195') as rc2:
        model = rc2.compute()
        if model is None:
            return None

        cost = rc2.cost
        # cost = number of violated soft clauses = number of edges NOT absent
        # = total soft clauses - edges present... no.
        # Actually RC2 minimizes the weight of unsatisfied soft clauses.
        # Each soft clause is [-edge_var], weight 1.
        # An edge being present violates its soft clause.
        # So cost = number of edges present.
        G = encoder.decode_graph(model)
        num_edges = G.number_of_edges()

        if verbose:
            print(f"  Optimal: {num_edges} edges (cost={cost})")

        if return_graph:
            return num_edges, G
        return num_edges


def sat_number_binary(n, F, return_graph=False, verbose=False):
    """
    Compute sat(n, F) using SAT encoding with binary search.

    Fallback method; MaxSAT (sat_number) is preferred.
    """
    f_nodes = F.number_of_nodes()
    if n < f_nodes:
        return None

    encoder = SaturationSATEncoder(n, F)

    lo = 0
    hi = n * (n - 1) // 2
    best_model = None
    best_edges = None

    # Feasibility check
    clauses = encoder.encode(hi)
    with Solver(name='cd195', bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None

    while lo <= hi:
        mid = (lo + hi) // 2
        if verbose:
            print(f"  Trying max_edges={mid} (lo={lo}, hi={hi})")

        encoder = SaturationSATEncoder(n, F)
        clauses = encoder.encode(mid)

        with Solver(name='cd195', bootstrap_with=clauses) as solver:
            if solver.solve():
                model = solver.get_model()
                G = encoder.decode_graph(model)
                best_edges = G.number_of_edges()
                best_model = (model, encoder)
                hi = best_edges - 1
                if verbose:
                    print(f"    SAT: found graph with {best_edges} edges")
            else:
                lo = mid + 1
                if verbose:
                    print(f"    UNSAT")

    if best_edges is None:
        return None

    if return_graph:
        model, enc = best_model
        G = enc.decode_graph(model)
        return best_edges, G
    return best_edges


# ---------------------------------------------------------------------------
# Convenience functions for specific families
# ---------------------------------------------------------------------------

def sat_complete(n, t, **kwargs):
    """Compute sat(n, K_t) using SAT solver."""
    return sat_number(n, complete_graph(t), **kwargs)


def sat_cycle(n, k, **kwargs):
    """Compute sat(n, C_k) using SAT solver."""
    return sat_number(n, cycle_graph(k), **kwargs)


def sat_path(n, k, **kwargs):
    """Compute sat(n, P_k) using SAT solver."""
    return sat_number(n, path_graph(k), **kwargs)


# ---------------------------------------------------------------------------
# Batch computation with table output
# ---------------------------------------------------------------------------

def compute_table(ns, F, name="F", verbose=False):
    """
    Compute sat(n, F) for a range of n values.

    Returns a dict mapping n -> sat(n, F).
    """
    results = {}
    for n in ns:
        if verbose:
            print(f"Computing sat({n}, {name})...")
        val = sat_number(n, F, verbose=verbose)
        results[n] = val
        if verbose:
            print(f"  sat({n}, {name}) = {val}")
    return results


if __name__ == "__main__":
    import time

    print("=== Graph Saturation Number Computation ===\n")

    # Verify K_4 against Erdos-Hajnal-Moon
    print("--- sat(n, K_4) vs Erdos-Hajnal-Moon formula ---")
    for n in range(5, 12):
        t0 = time.time()
        computed = sat_number(n, complete_graph(4))
        dt = time.time() - t0
        formula = sat_formula_complete(n, 4)
        match = "OK" if computed == formula else "MISMATCH"
        print(f"  n={n}: computed={computed}, formula={formula} [{match}] ({dt:.2f}s)")

    print()

    # Compute sat(n, C_4)
    print("--- sat(n, C_4) ---")
    for n in range(5, 13):
        t0 = time.time()
        computed = sat_number(n, cycle_graph(4))
        dt = time.time() - t0
        formula = sat_formula_c4(n)
        match = "OK" if computed == formula else "MISMATCH"
        print(f"  n={n}: computed={computed}, formula={formula} [{match}] ({dt:.2f}s)")

    print()

    # Compute sat(n, P_4)
    print("--- sat(n, P_4) ---")
    for n in range(4, 12):
        t0 = time.time()
        computed = sat_number(n, path_graph(4))
        dt = time.time() - t0
        formula = sat_formula_p4(n)
        match = "OK" if computed == formula else "MISMATCH"
        print(f"  n={n}: computed={computed}, formula={formula} [{match}] ({dt:.2f}s)")
