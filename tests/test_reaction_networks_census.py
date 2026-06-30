"""Tests for reaction_networks_census.py.

Validation strategy (verification-always, exact arithmetic):

  1. ``integer_rank`` matches ``sympy.Matrix.rank`` bit-for-bit on a large
     sweep of stoichiometry matrices (the fast path must equal the slow path).
  2. The fast invariants (``deficiency_fast``, ``is_weakly_reversible_fast``)
     agree with the reference implementations in ``reaction_networks`` on every
     enumerated small network.
  3. The census reproduces, exactly, the hand-verified cases in
     ``discoveries/biology/reaction-networks/findings.md`` and the per-r
     regression values frozen in ``test_reaction_networks.py``.
  4. The orbit counts produced by the census agree with an *independent*
     Burnside's-lemma count.
  5. Structural invariants: ``delta >= 0`` always; ``delta == 0`` iff
     ``s == n - ell``; weakly reversible deficiency-zero networks satisfy the
     hypotheses of Feinberg's Deficiency Zero Theorem.
  6. Frozen complete-census values (k=1 exhaustive, k=2 complete 2^20).
"""

from __future__ import annotations

import pytest
from sympy import Matrix

from src.reaction_networks import (
    Network,
    deficiency,
    enumerate_networks,
    enumerate_reactions,
    is_weakly_reversible,
    linkage_classes,
    stoichiometry_rank,
    tabulate,
)
from src.reaction_networks_census import (
    CensusResult,
    census,
    deficiency_fast,
    integer_rank,
    is_weakly_reversible_fast,
    labelled_total,
    orbit_count_burnside,
)


# ---------------------------------------------------------------------------
# 1. integer_rank == sympy rank
# ---------------------------------------------------------------------------


def test_integer_rank_matches_sympy_small_cases():
    cases = [
        [[1, 0], [0, 1]],            # rank 2
        [[1, -1], [1, -1]],          # rank 1 (collinear, Feinberg p.6)
        [[1, 0, 0], [0, 0, 0]],      # rank 1
        [[2, -2, 0], [-1, 1, 0]],    # rank 1
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # rank 3
        [[1, 2, 3], [2, 4, 6], [1, 1, 1]],  # rank 2
        [],                          # rank 0
    ]
    for rows in cases:
        expected = int(Matrix(rows).rank()) if rows else 0
        assert integer_rank([list(r) for r in rows]) == expected, rows


def test_integer_rank_matches_sympy_on_enumerated_networks():
    """Bit-for-bit agreement on every k=3, r<=3 network's stoichiometry matrix."""
    checked = 0
    for net in enumerate_networks(k=3, r=3, max_stoich=2):
        rows = [[yp[i] - y[i] for i in range(3)] for (y, yp) in net.reactions]
        fast = integer_rank([list(r) for r in rows])
        slow = stoichiometry_rank(net)
        assert fast == slow, (net.reactions, fast, slow)
        checked += 1
    assert checked > 1000


# ---------------------------------------------------------------------------
# 2. fast invariants == reference invariants
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k,r", [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2)])
def test_fast_invariants_match_reference(k, r):
    for net in enumerate_networks(k=k, r=r, max_stoich=2):
        rx = net.reactions
        assert deficiency_fast(rx, k) == deficiency(net)
        assert is_weakly_reversible_fast(rx) == is_weakly_reversible(net)


# ---------------------------------------------------------------------------
# 3a. hand-verified cases from findings.md (n, ell, s, delta)
# ---------------------------------------------------------------------------


def _nlsd(reactions: tuple, k: int) -> tuple[int, int, int, int]:
    net = Network(k=k, reactions=tuple(sorted(reactions)))
    n = len(net.complexes())
    ell = len(linkage_classes(net))
    s = stoichiometry_rank(net)
    d = deficiency_fast(reactions, k)
    return (n, ell, s, d)


def test_findings_hand_verified_cases():
    A1, B1 = (1, 0), (0, 1)
    # A -> B : (2,1,1,0)
    assert _nlsd(((A1, B1),), 2) == (2, 1, 1, 0)
    # A <-> B : (2,1,1,0), WR
    assert _nlsd(((A1, B1), (B1, A1)), 2) == (2, 1, 1, 0)
    assert is_weakly_reversible_fast(((A1, B1), (B1, A1)))
    # A + B -> C (k=3): (2,1,1,0)
    AB, C = (1, 1, 0), (0, 0, 1)
    assert _nlsd(((AB, C),), 3) == (2, 1, 1, 0)
    # A -> 2B, 2B -> A+C : (3,1,2,0)
    A3, B2, AC = (1, 0, 0), (0, 2, 0), (1, 0, 1)
    assert _nlsd(((A3, B2), (B2, AC)), 3) == (3, 1, 2, 0)
    # A+B -> 2A, B -> A : (4,2,1,1)  [Feinberg p.6]
    ABx, AA, Bx, Ax = (1, 1), (2, 0), (0, 1), (1, 0)
    assert _nlsd(((ABx, AA), (Bx, Ax)), 2) == (4, 2, 1, 1)
    # A -> B -> C -> A : (3,1,2,0), WR
    Av, Bv, Cv = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    cyc = ((Av, Bv), (Bv, Cv), (Cv, Av))
    assert _nlsd(cyc, 3) == (3, 1, 2, 0)
    assert is_weakly_reversible_fast(cyc)


# ---------------------------------------------------------------------------
# 3b. census reproduces per-r regression values (frozen by original module)
# ---------------------------------------------------------------------------


def test_census_by_r_matches_tabulate_k2():
    res = census(k=2, max_stoich=2, r_max=3)
    for r in (1, 2, 3):
        t = tabulate(k=2, r=r, max_stoich=2)
        c = res.by_r[r]
        assert c["total"] == t["total"]
        assert c["by_deficiency"] == t["by_deficiency"]
        assert c["wr"] == t["weakly_reversible_count"]


def test_census_by_r_matches_findings_k3_small():
    """findings.md k=3 table: r=1->5, r=2->350, with deficiency splits."""
    res = census(k=3, max_stoich=2, r_max=2)
    assert res.by_r[1]["total"] == 5
    assert res.by_r[1]["by_deficiency"] == {0: 5}
    assert res.by_r[2]["total"] == 350
    assert res.by_r[2]["by_deficiency"] == {0: 339, 1: 11}


def test_census_per_r_full_k2_table_from_findings():
    """The complete k=2 per-r table (r=1..10) from findings.md, exactly.

    Slow-ish (covers r up to 6 brute combinations); the r=7..10 rows are
    checked in the slow complete-census test below.
    """
    expected = {
        1: {"total": 8, "by_deficiency": {0: 8}, "wr": 0},
        2: {"total": 99, "by_deficiency": {0: 82, 1: 17}, "wr": 5},
        3: {"total": 570, "by_deficiency": {0: 110, 1: 448, 2: 12}, "wr": 10},
        4: {"total": 2445, "by_deficiency": {0: 75, 1: 1333, 2: 1037}, "wr": 71},
        5: {"total": 7752, "by_deficiency": {0: 27, 1: 2148, 2: 5577}, "wr": 262},
        6: {"total": 19440, "by_deficiency": {0: 5, 1: 2434, 2: 17001}, "wr": 1210},
    }
    res = census(k=2, max_stoich=2, r_max=6)
    for r, exp in expected.items():
        c = res.by_r[r]
        assert c["total"] == exp["total"], (r, c["total"])
        assert c["by_deficiency"] == exp["by_deficiency"], (r, c["by_deficiency"])
        assert c["wr"] == exp["wr"], (r, c["wr"])


# ---------------------------------------------------------------------------
# 4. Burnside cross-check on orbit counts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [1, 2, 3])
def test_burnside_matches_census_orbits(k):
    res = census(k=k, max_stoich=2, r_max=3)
    for r in (1, 2, 3):
        burn = orbit_count_burnside(k, 2, r)
        orb = res.by_r.get(r, {}).get("total", 0)
        assert burn == orb, (k, r, burn, orb)


def test_burnside_labelled_relation_k2_r2():
    """Sanity: orbit count <= labelled count, and Burnside is an integer avg."""
    lab = labelled_total(2, 2, 2)
    orb = orbit_count_burnside(2, 2, 2)
    assert orb <= lab
    assert orb == 99  # matches census


# ---------------------------------------------------------------------------
# 5. Structural invariants (Feinberg)
# ---------------------------------------------------------------------------


def test_deficiency_nonnegative_everywhere_k2_complete_small():
    """delta >= 0 for every network in the k=2 census up to r=4."""
    res = census(k=2, max_stoich=2, r_max=4)
    assert min(res.by_deficiency) >= 0


def test_deficiency_zero_iff_rank_equals_n_minus_ell():
    """delta == 0  <=>  s == n - ell, checked on every k=3, r<=3 network.

    This is the structural identity underlying the Deficiency Zero Theorem:
    the stoichiometric subspace has the maximal dimension n - ell exactly when
    the deficiency vanishes.
    """
    for net in enumerate_networks(k=3, r=3, max_stoich=2):
        n = len(net.complexes())
        ell = len(linkage_classes(net))
        s = stoichiometry_rank(net)
        d = deficiency_fast(net.reactions, 3)
        assert (d == 0) == (s == n - ell)


def test_weakly_reversible_deficiency_zero_satisfies_dzt_hypotheses():
    """Every (WR, delta=0) network in a small census meets DZT preconditions.

    Feinberg's Deficiency Zero Theorem requires: (i) delta = 0, (ii) weak
    reversibility.  Under those hypotheses the network must also have
    s = n - ell (structural) and every linkage class strongly connected
    (definition of WR).  We confirm both hold for the entire (WR, delta=0)
    subpopulation, and that the subpopulation is non-empty.
    """
    found = 0
    for net in enumerate_networks(k=2, r=2, max_stoich=2):
        rx = net.reactions
        if deficiency_fast(rx, 2) == 0 and is_weakly_reversible_fast(rx):
            found += 1
            n = len(net.complexes())
            ell = len(linkage_classes(net))
            s = stoichiometry_rank(net)
            assert s == n - ell
            # WR <=> reference WR (every linkage class strongly connected)
            assert is_weakly_reversible(net)
    assert found == 5  # findings.md: 5 such networks at k=2, r=2


def test_deficiency_zero_examples_are_deficiency_zero():
    """Spot-check named delta=0 networks via the fast path."""
    A, B = (1, 0), (0, 1)
    assert deficiency_fast(((A, B),), 2) == 0
    assert deficiency_fast(((A, B), (B, A)), 2) == 0


# ---------------------------------------------------------------------------
# 6. Complete census: frozen exact values
# ---------------------------------------------------------------------------


def test_complete_census_k1():
    """k=1, max_stoich=2: only complexes A=(1,), 2A=(2,).

    Networks (all species used, no zero complex): {A->2A}, {2A->A},
    {A<->2A}.  All deficiency 0; the reversible one is weakly reversible.
    S_1 is trivial so A->2A and 2A->A are distinct.
    """
    res = census(k=1, max_stoich=2)
    assert res.complete
    assert res.total == 3
    assert res.by_deficiency == {0: 3}
    assert res.weakly_reversible == 1
    assert res.wr_and_def0 == 1
    assert res.by_r[1]["total"] == 2
    assert res.by_r[2]["total"] == 1


@pytest.mark.slow
def test_complete_census_k2_full_2pow20():
    """Complete k=2 census over all 2^20 reaction subsets (frozen values).

    Total networks (any r, up to S_2) and the full deficiency distribution,
    weak reversibility, and WR&delta=0 counts.  Also re-confirms the entire
    per-r table r=1..20 (extends the original findings.md table, which
    stopped at r=10).
    """
    res = census(k=2, max_stoich=2, r_max=None, workers=None)
    assert res.complete
    assert res.total == 524796
    assert res.by_deficiency == {0: 307, 1: 10397, 2: 514092}
    assert res.weakly_reversible == 286982
    assert res.wr_and_def0 == 93

    # full per-r table (r : (total, by_deficiency, wr))
    per_r = {
        1: (8, {0: 8}, 0),
        2: (99, {0: 82, 1: 17}, 5),
        3: (570, {0: 110, 1: 448, 2: 12}, 10),
        4: (2445, {0: 75, 1: 1333, 2: 1037}, 71),
        5: (7752, {0: 27, 1: 2148, 2: 5577}, 262),
        6: (19440, {0: 5, 1: 2434, 2: 17001}, 1210),
        7: (38760, {1: 2016, 2: 36744}, 4480),
        8: (63090, {1: 1250, 2: 61840}, 14510),
        9: (83980, {1: 550, 2: 83430}, 33170),
        10: (92504, {1: 168, 2: 92336}, 53000),
        11: (83980, {1: 30, 2: 83950}, 61320),
        12: (63090, {1: 3, 2: 63087}, 53482),
        13: (38760, {2: 38760}, 35850),
        14: (19440, {2: 19440}, 18820),
        15: (7752, {2: 7752}, 7672),
        16: (2445, {2: 2445}, 2439),
        17: (570, {2: 570}, 570),
        18: (100, {2: 100}, 100),
        19: (10, {2: 10}, 10),
        20: (1, {2: 1}, 1),
    }
    for r, (tot, bd, wr) in per_r.items():
        c = res.by_r[r]
        assert c["total"] == tot, (r, c["total"])
        assert c["by_deficiency"] == bd, (r, c["by_deficiency"])
        assert c["wr"] == wr, (r, c["wr"])

    # cross-check: sum of per-r totals == grand total
    assert sum(res.by_r[r]["total"] for r in res.by_r) == res.total


# ---------------------------------------------------------------------------
# CensusResult helpers
# ---------------------------------------------------------------------------


def test_census_result_sequence_helpers():
    res = census(k=2, max_stoich=2, r_max=3)
    assert isinstance(res, CensusResult)
    seq = res.deficiency_sequence()
    # cumulative over r<=3: delta 0,1,2 present
    assert len(seq) == 3
    assert all(x >= 0 for x in seq)
    rseq = res.r_total_sequence()
    assert rseq[:3] == [8, 99, 570]
