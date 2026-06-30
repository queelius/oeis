#!/usr/bin/env python3
"""Reproduce and verify A395521.

a(n) = number of non-isomorphic abelian groups appearing as the sandpile
group K(G) (critical group / Jacobian) of a connected graph G on n vertices.

K(G) = Z^(n-1) / image(reduced Laplacian); its invariant factors come from
the Smith normal form of the reduced Laplacian.  |K(G)| equals the number
of spanning trees of G (matrix-tree theorem).

Generator: src/sandpile_groups/build_catalog.py, driven by src/sandpile.py
  - enumerate_connected_graphs(n): every isomorphism class of connected
        simple graph on n vertices (pynauty canonical certificates).
  - sandpile_group(adj) -> invariant factors via Smith normal form.
  - group_signature(...) -> the abstract-group key 'Z/d1 x Z/d2 x ...';
        distinct signatures over all G give a(n).

This reproduce.py recomputes the ENTIRE published DATA (n = 2..8) from
scratch (~25 s total) and additionally cross-checks the connected-graph
count at each n against OEIS A001349.

Run as:  python3 reproduce.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 "..", "..")))

from src.sandpile import (  # noqa: E402
    adjacency_matrix,
    enumerate_connected_graphs,
    group_signature,
    sandpile_group,
)

ID = "A395521"
OFFSET = 2
# Published DATA a(n) for n = 2..8.
DATA = [1, 2, 5, 16, 71, 461, 4855]
# OEIS A001349: number of connected simple graphs on n nodes (cross-check).
A001349 = {2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117}


def distinct_sandpile_groups(n: int) -> tuple[int, int]:
    """Return (#connected graphs on n vertices, #distinct sandpile groups)."""
    signatures: set[str] = set()
    graphs = 0
    for (_, edges) in enumerate_connected_graphs(n):
        graphs += 1
        adj = adjacency_matrix(n, edges)
        K = sandpile_group(adj, root=0)
        signatures.add(group_signature(K["invariant_factors"]))
    return graphs, len(signatures)


def main() -> None:
    n_terms = len(DATA)
    for i, expected in enumerate(DATA):
        n = OFFSET + i
        graphs, distinct = distinct_sandpile_groups(n)
        assert graphs == A001349[n], \
            f"connected graphs on {n} = {graphs} != A001349 {A001349[n]}"
        assert distinct == expected, \
            f"distinct sandpile groups on {n} = {distinct} != a({n})={expected}"

    print(f"OK {ID}: reproduced {n_terms} terms (n=2..{OFFSET + n_terms - 1}) "
          f"from scratch via connected-graph enumeration + Smith normal form; "
          f"connected-graph counts match A001349")


if __name__ == "__main__":
    main()
