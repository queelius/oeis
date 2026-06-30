#!/usr/bin/env python3
"""
Fast backtracking solver for Schur numbers S(G, k) in finite abelian groups.

This module is the value-generating code for the k=3 cyclic Schur sequence
S(Z/nZ, 3).  The brute-force `schur_number(..., k=3)` in `schur_groups`
enumerates all k^|subset| colorings and becomes intractable past n ~ 11;
`backtrack_schur` below uses backtracking with forward checking on the fixed
Schur-triple constraint graph and reaches n = 15 essentially instantly.

Both `_schur_triples` and `backtrack_schur` are copied verbatim (modulo this
docstring) from the source repo
    github.com/queelius/computational-explorations  ->  src/schur_algorithms.py
The original `schur_algorithms.py` imports `pysat` at module top for its
optional SAT-based solvers; those are NOT needed for the backtracking solver
used here, so this companion module deliberately omits the SAT machinery and
has no third-party dependencies (stdlib only, plus the sibling
`schur_groups`).

S(G, k) = the maximum size of a subset A of G that can be partitioned into k
sum-free sets (equivalently, k-colored with every color class sum-free).
"""

from itertools import combinations
from typing import Dict, List, Set, Tuple

try:  # when imported with the oeis repo root on sys.path (from src.<module>)
    from src.schur_groups import group_add, group_elements, is_sum_free
except ImportError:  # when src/ itself is on sys.path (from <module>)
    from schur_groups import group_add, group_elements, is_sum_free


def _schur_triples(elements: List[Tuple[int, ...]],
                   orders: Tuple[int, ...]) -> List[Tuple[int, int, int]]:
    """Return all (i, j, k) index-triples with elements[i]+elements[j]=elements[k].

    Only includes triples where i and j are both in the candidate subset
    (which will be constrained later). Precomputing this is the key to
    making SAT / backtracking fast: the constraint graph is fixed.
    """
    elem_to_idx = {e: i for i, e in enumerate(elements)}
    triples = []
    n = len(elements)
    for i in range(n):
        for j in range(i, n):  # j >= i avoids duplicate (i,j)/(j,i)
            s = group_add(elements[i], elements[j], orders)
            if s in elem_to_idx:
                k = elem_to_idx[s]
                triples.append((i, j, k))
                if i != j:
                    triples.append((j, i, k))
    return triples


def backtrack_schur(orders: Tuple[int, ...], k: int) -> int:
    """Compute S(G,k) via backtracking with constraint propagation.

    For each candidate subset (tried largest-first), attempts to
    k-color it using backtracking with forward checking:
    when element i is assigned color c, all Schur-constrained elements
    have color c removed from their domain.
    """
    elements = group_elements(orders)
    n = len(elements)
    if n > 30:
        return -1

    triples = _schur_triples(elements, orders)

    # Build adjacency: for each element, ALL triples it participates in,
    # regardless of position. Each entry stores the triple and the
    # OTHER two element indices.
    # Constraint: for triple (a, b, s) with a+b=s, all three in the
    # same color is forbidden. So when elem gets color c, for every
    # triple containing elem, if the other two elements both have color c,
    # that's a conflict.
    # Forward check: when elem gets color c, for each triple (a,b,s)
    # containing elem, check the other two elements (p, q).
    # If one of (p,q) already has color c, the other cannot have color c.
    elem_triples_full: Dict[int, List[Tuple[int, int]]] = {i: [] for i in range(n)}
    seen = set()
    for (a, b, s) in triples:
        triple_key = tuple(sorted([a, b, s]))
        if triple_key in seen:
            continue
        seen.add(triple_key)
        elem_triples_full[a].append((b, s))
        elem_triples_full[b].append((a, s))
        elem_triples_full[s].append((a, b))

    def _backtrack(assignment: List[int], domains: List[Set[int]],
                   idx: int, selected: List[int]) -> bool:
        """Try to color selected[idx..] given current assignment and domains."""
        if idx == len(selected):
            return True
        elem = selected[idx]
        for c in sorted(domains[elem]):
            # Forward check: would assigning elem=c empty any future domain?
            pruned: List[Tuple[int, int]] = []
            feasible = True
            assignment[elem] = c

            for (p, q) in elem_triples_full[elem]:
                # elem, p, q form a forbidden triple in color c.
                # If p has color c, q cannot have color c.
                # If q has color c, p cannot have color c.
                if p not in _sel_set or q not in _sel_set:
                    continue
                if assignment[p] == c and assignment[q] == -1:
                    if c in domains[q]:
                        domains[q].discard(c)
                        pruned.append((q, c))
                        if not domains[q]:
                            feasible = False
                            break
                elif assignment[q] == c and assignment[p] == -1:
                    if c in domains[p]:
                        domains[p].discard(c)
                        pruned.append((p, c))
                        if not domains[p]:
                            feasible = False
                            break
                elif assignment[p] == c and assignment[q] == c:
                    # Conflict: all three have color c
                    feasible = False
                    break

            if feasible and _backtrack(assignment, domains, idx + 1, selected):
                return True

            # Undo
            assignment[elem] = -1
            for (e, cv) in pruned:
                domains[e].add(cv)

        return False

    # Try subsets from largest down
    best = 0
    zero_elem = elements.index(tuple(0 for _ in orders))
    for size in range(n, 0, -1):
        if size <= best:
            break
        for combo in combinations(range(n), size):
            selected = list(combo)
            _sel_set = set(selected)

            # Zero can never be in a sum-free set (0+0=0); skip
            if zero_elem in _sel_set:
                continue

            assignment = [-1] * n
            domains = [set(range(k)) if i in _sel_set else set()
                       for i in range(n)]

            if _backtrack(assignment, domains, 0, selected):
                return size
    return best
