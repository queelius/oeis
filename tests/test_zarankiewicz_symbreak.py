"""Tests for double-lex symmetry-broken Zarankiewicz encoding.

Two layers:
1. The lex_leq CNF is logically correct (brute-forced on small vectors).
2. The symmetry break is SOUND: it reproduces known exact z values
   (SAT at the value, UNSAT at value+1) -- same answers as the plain
   encoding, which is the whole point.
"""

from __future__ import annotations

import itertools

import pytest

from src.zarankiewicz import is_kst_free, count_ones
from src.zarankiewicz_symbreak import (
    build_encoding,
    decode_matrix,
    lex_leq_clauses,
    double_lex_clauses,
)

pysat = pytest.importorskip("pysat")
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType


# ---- lex_leq correctness (brute force) ----

def _lex_leq_truth(a_bits, b_bits):
    return list(a_bits) <= list(b_bits)


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_lex_leq_matches_truth_table(k):
    # Variables 1..k = a, k+1..2k = b. Enumerate all assignments, check the
    # CNF is satisfiable-consistent with a <=_lex b exactly.
    a = list(range(1, k + 1))
    b = list(range(k + 1, 2 * k + 1))
    clauses, top = lex_leq_clauses(a, b, 2 * k + 1)
    for a_bits in itertools.product([0, 1], repeat=k):
        for b_bits in itertools.product([0, 1], repeat=k):
            # Build assignment for a, b; let solver choose aux vars.
            assumptions = []
            for idx, bit in enumerate(a_bits):
                assumptions.append(a[idx] if bit else -a[idx])
            for idx, bit in enumerate(b_bits):
                assumptions.append(b[idx] if bit else -b[idx])
            with Solver(name="g3", bootstrap_with=clauses) as s:
                sat = s.solve(assumptions=assumptions)
            assert sat == _lex_leq_truth(a_bits, b_bits), (a_bits, b_bits, sat)


def test_double_lex_clause_vars_disjoint_from_matrix():
    # Aux vars must start above m*n.
    m, n = 4, 5
    clauses, top = double_lex_clauses(m, n, m * n + 1)
    aux_used = {abs(l) for cl in clauses for l in cl if abs(l) > m * n}
    assert all(v > m * n for v in aux_used)


# ---- soundness: reproduce known exact z(m,n;3,4) values ----

# Known exact values (from src/zarankiewicz KNOWN tables).
KNOWN = {
    (9, 9, 3, 4): 56,
    (9, 10, 3, 4): 60,
    (9, 11, 3, 4): 66,
    (10, 10, 3, 4): 66,
}


def _decide(m, n, s, t, weight, symmetry, budget_solver="g3"):
    num_vars, clauses, top_id = build_encoding(m, n, s, t, symmetry=symmetry)
    card = CardEnc.atleast(lits=list(range(1, num_vars + 1)), bound=weight,
                           top_id=top_id, encoding=EncType.seqcounter)
    with Solver(name=budget_solver, bootstrap_with=clauses) as s_:
        for cl in card.clauses:
            s_.add_clause(cl)
        sat = s_.solve()
        model = s_.get_model() if sat else None
    return sat, model


@pytest.mark.slow
@pytest.mark.parametrize("m,n,s,t,val", [(k[0], k[1], k[2], k[3], v) for k, v in KNOWN.items()])
def test_symbreak_reproduces_known_sat_at_value(m, n, s, t, val):
    # SAT at the exact value: a doubly-lex witness exists.
    sat, model = _decide(m, n, s, t, val, symmetry=True)
    assert sat, f"z({m},{n};{s},{t}) should be SAT at weight {val} under symmetry"
    M = decode_matrix(model, m, n)
    assert is_kst_free(M, s, t)
    assert count_ones(M) >= val


@pytest.mark.slow
@pytest.mark.parametrize("m,n,s,t,val", [(k[0], k[1], k[2], k[3], v) for k, v in KNOWN.items()])
def test_symbreak_reproduces_known_unsat_above(m, n, s, t, val):
    # UNSAT at value+1: this is the upper-bound proof the symmetry break enables.
    sat, _ = _decide(m, n, s, t, val + 1, symmetry=True)
    assert not sat, f"z({m},{n};{s},{t}) should be UNSAT at weight {val+1} under symmetry"


def test_symbreak_small_consistency_with_plain():
    # On a small case, symmetry-broken and plain must agree on SAT/UNSAT.
    m, n, s, t, w = 9, 9, 3, 4, 56
    sat_sym, _ = _decide(m, n, s, t, w, symmetry=True)
    sat_plain, _ = _decide(m, n, s, t, w, symmetry=False)
    assert sat_sym == sat_plain == True
