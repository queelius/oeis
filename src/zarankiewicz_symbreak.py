"""Zarankiewicz SAT with double-lex symmetry breaking.

The plain z(m,n;s,t) encoding (src/zarankiewicz.py) forbids every s x t
all-ones submatrix and adds a cardinality bound. It has NO symmetry
breaking, so the solver must rule out all m! row-permutations and n!
column-permutations of every candidate -- which makes the UNSAT side
(proving an upper bound z < w) intractable at the phase transition.

This module adds DOUBLE-LEX symmetry breaking: rows constrained to be in
non-decreasing lexicographic order, and columns likewise. This is the
standard sound symmetry break for matrix models with independent row and
column symmetry (Flener, Frisch, Hnich, Kiziltan, Miguel, Pearson, Walsh,
"Breaking row and column symmetries in matrix models", CP 2002): every
orbit under row x column permutations contains at least one doubly-lex
matrix, so satisfiability is preserved while the search space collapses.

The lex constraint a <=_lex b between two equal-length bit vectors uses
the standard prefix-equality auxiliary encoding (Knuth TAOCP 7.2.2.2;
Codish et al.).

Soundness: a witness exists under double-lex iff one exists without it
(row/col permutations preserve both weight and K_{s,t}-freeness). So:
  * SAT at w   <=>  z(m,n;s,t) >= w   (same as plain encoding)
  * UNSAT at w <=>  z(m,n;s,t) <  w   (same, but reachable)
"""

from __future__ import annotations

from itertools import combinations


def _var(i: int, j: int, n: int) -> int:
    """Matrix entry (i, j) -> SAT variable (1-indexed), row-major."""
    return i * n + j + 1


def build_forbidden_clauses(m: int, n: int, s: int, t: int) -> tuple[int, list[list[int]]]:
    """Clauses forbidding every s x t all-ones submatrix. (m*n vars.)"""
    clauses = []
    for row_subset in combinations(range(m), s):
        for col_subset in combinations(range(n), t):
            clauses.append([-_var(r, c, n) for r in row_subset for c in col_subset])
    return m * n, clauses


def lex_leq_clauses(a: list[int], b: list[int], next_var: int) -> tuple[list[list[int]], int]:
    """CNF enforcing a <=_lex b for equal-length literal lists a, b.

    Uses a running 'prefix-equal' literal. Returns (clauses, next_var).
    """
    assert len(a) == len(b)
    clauses: list[list[int]] = []
    k = len(a)
    prefix: int | None = None  # p_j; None denotes constant-true (empty prefix equal)

    for j in range(k):
        aj, bj = a[j], b[j]
        # If equal so far, require a_j <= b_j, i.e. (¬prefix ∨ ¬a_j ∨ b_j).
        if prefix is None:
            clauses.append([-aj, bj])
        else:
            clauses.append([-prefix, -aj, bj])

        if j == k - 1:
            break  # no prefix needed past the last position

        # eq_j <-> (a_j == b_j)
        eqj = next_var
        next_var += 1
        clauses.append([-eqj, -aj, bj])   # eq -> (a -> b)
        clauses.append([-eqj, -bj, aj])   # eq -> (b -> a)
        clauses.append([aj, bj, eqj])     # (¬a ∧ ¬b) -> eq
        clauses.append([-aj, -bj, eqj])   # (a ∧ b) -> eq

        # p_{j+1} <-> prefix ∧ eq_j
        if prefix is None:
            prefix = eqj  # p_{j+1} == eq_j  (since p_j is true)
        else:
            pj1 = next_var
            next_var += 1
            clauses.append([-pj1, prefix])         # p_{j+1} -> prefix
            clauses.append([-pj1, eqj])            # p_{j+1} -> eq_j
            clauses.append([-prefix, -eqj, pj1])   # prefix ∧ eq_j -> p_{j+1}
            prefix = pj1

    return clauses, next_var


def double_lex_clauses(m: int, n: int, next_var: int) -> tuple[list[list[int]], int]:
    """Row-lex (row i <= row i+1) and column-lex (col j <= col j+1) clauses."""
    clauses: list[list[int]] = []
    # Rows non-decreasing in lex order.
    for i in range(m - 1):
        a = [_var(i, j, n) for j in range(n)]
        b = [_var(i + 1, j, n) for j in range(n)]
        cl, next_var = lex_leq_clauses(a, b, next_var)
        clauses.extend(cl)
    # Columns non-decreasing in lex order.
    for j in range(n - 1):
        a = [_var(i, j, n) for i in range(m)]
        b = [_var(i, j + 1, n) for i in range(m)]
        cl, next_var = lex_leq_clauses(a, b, next_var)
        clauses.extend(cl)
    return clauses, next_var


def build_encoding(m: int, n: int, s: int, t: int, symmetry: bool = True):
    """Full clause list: forbidden + (optional) double-lex. Returns (num_vars, clauses, top_id)."""
    num_vars, clauses = build_forbidden_clauses(m, n, s, t)
    top_id = num_vars
    if symmetry:
        sym_clauses, top_id = double_lex_clauses(m, n, num_vars + 1)
        clauses = clauses + sym_clauses
        top_id = max(top_id - 1, num_vars)
    return num_vars, clauses, top_id


def decode_matrix(model: list[int], m: int, n: int) -> list[list[int]]:
    true_set = {v for v in model if v > 0}
    return [[1 if _var(i, j, n) in true_set else 0 for j in range(n)] for i in range(m)]
