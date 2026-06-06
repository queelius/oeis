"""Tests for reaction_networks.py."""

from __future__ import annotations

import pytest

from src.reaction_networks import (
    Network,
    canonical_form,
    deficiency,
    enumerate_complexes,
    enumerate_networks,
    enumerate_reactions,
    is_deficiency_zero,
    is_weakly_reversible,
    linkage_classes,
    stoichiometry_rank,
    tabulate,
)


# ---------------------------------------------------------------------------
# Basic invariants
# ---------------------------------------------------------------------------


def test_a_to_b_is_deficiency_zero():
    """A -> B has 2 complexes, 1 linkage class, rank 1, delta = 0."""
    A, B = (1, 0), (0, 1)
    net = Network(k=2, reactions=((A, B),))
    assert len(net.complexes()) == 2
    assert len(linkage_classes(net)) == 1
    assert stoichiometry_rank(net) == 1
    assert deficiency(net) == 0
    assert is_deficiency_zero(net)
    assert not is_weakly_reversible(net)  # one-way reaction


def test_a_to_b_b_to_a_is_weakly_reversible():
    """A <-> B is weakly reversible with delta = 0."""
    A, B = (1, 0), (0, 1)
    net = Network(k=2, reactions=((A, B), (B, A)))
    assert deficiency(net) == 0
    assert is_weakly_reversible(net)


def test_apb_to_c_three_species():
    """A + B -> C: 2 complexes (A+B, C), 1 linkage class, rank 1, delta=0."""
    AB = (1, 1, 0)
    C = (0, 0, 1)
    net = Network(k=3, reactions=((AB, C),))
    assert len(net.complexes()) == 2
    assert deficiency(net) == 0


def test_feinberg_classic_deficiency_one():
    """A+B -> 2A, B -> A: 4 complexes, 2 linkage classes, rank 1, delta=1.

    From Feinberg's textbook page 6.  Stoichiometric vectors are
    2A - (A+B) = (1,-1) and A - B = (1,-1), which span a 1-d subspace.
    """
    AB = (1, 1)
    AA = (2, 0)
    B = (0, 1)
    A = (1, 0)
    net = Network(k=2, reactions=((AB, AA), (B, A)))
    assert len(net.complexes()) == 4
    assert len(linkage_classes(net)) == 2
    assert stoichiometry_rank(net) == 1
    assert deficiency(net) == 1


def test_horn_jackson_deficiency_one():
    """A -> 2B, 2B -> A+C: 3 complexes, 1 linkage class, rank 2, delta=0.

    The two stoichiometric vectors 2B-A=(-1,2,0) and (A+C)-2B=(1,-2,1) are
    linearly independent, so s=2 and delta = 3 - 1 - 2 = 0.
    """
    A = (1, 0, 0)
    B2 = (0, 2, 0)
    AC = (1, 0, 1)
    net = Network(k=3, reactions=((A, B2), (B2, AC)))
    assert deficiency(net) == 0


def test_three_cycle_is_weakly_reversible():
    """A -> B -> C -> A is weakly reversible (single 3-cycle SCC)."""
    A = (1, 0, 0)
    B = (0, 1, 0)
    C = (0, 0, 1)
    net = Network(k=3, reactions=((A, B), (B, C), (C, A)))
    assert is_weakly_reversible(net)
    # 3 complexes, 1 linkage class, rank 2, delta = 0.
    assert deficiency(net) == 0


def test_chain_not_weakly_reversible():
    """A -> B -> C is *not* weakly reversible (no cycle)."""
    A = (1, 0, 0)
    B = (0, 1, 0)
    C = (0, 0, 1)
    net = Network(k=3, reactions=((A, B), (B, C)))
    assert not is_weakly_reversible(net)


def test_deficiency_nonnegative_on_random_small_networks():
    """Feinberg's theorem: delta >= 0 always.  Sweep all small networks."""
    for net in enumerate_networks(k=2, r=3, max_stoich=2):
        d = deficiency(net)
        assert d >= 0


# ---------------------------------------------------------------------------
# Canonicalization
# ---------------------------------------------------------------------------


def test_canonical_form_is_idempotent():
    A, B = (1, 0), (0, 1)
    net = Network(k=2, reactions=((A, B), (B, A)))
    c1 = canonical_form(net)
    c2 = canonical_form(c1)
    assert c1 == c2


def test_canonical_form_identifies_relabelings():
    """Two networks differing by species relabeling have the same canonical form."""
    A, B = (1, 0), (0, 1)
    net1 = Network(k=2, reactions=((A, B),))
    # Same network with species labels swapped: B -> A.
    net2 = Network(k=2, reactions=((B, A),))
    assert canonical_form(net1) == canonical_form(net2)


# ---------------------------------------------------------------------------
# Enumeration
# ---------------------------------------------------------------------------


def test_enumerate_complexes_k2_max_stoich_2():
    """6 complexes on 2 species with total <= 2: 0, A, B, 2A, A+B, 2B."""
    cs = enumerate_complexes(k=2, max_stoich=2)
    expected = {(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)}
    assert set(cs) == expected


def test_enumerate_reactions_k2_max_stoich_1_no_zero():
    """Only 2 reactions: A -> B and B -> A."""
    rs = enumerate_reactions(k=2, max_stoich=1, allow_zero_complex=False)
    assert set(rs) == {((1, 0), (0, 1)), ((0, 1), (1, 0))}


def test_enumerate_networks_k2_r1_count_is_eight():
    """Hand-counted: 8 inequivalent r=1 networks on k=2 with max_stoich=2."""
    nets = list(enumerate_networks(k=2, r=1, max_stoich=2))
    assert len(nets) == 8


def test_enumerate_networks_no_duplicates():
    """No two networks in the enumeration share a canonical form."""
    nets = list(enumerate_networks(k=2, r=2, max_stoich=2))
    sigs = {tuple(sorted(n.reactions)) for n in nets}
    assert len(sigs) == len(nets)


# ---------------------------------------------------------------------------
# Tabulation regression tests (frozen by initial sweep)
# ---------------------------------------------------------------------------


def test_tabulate_k2_r1():
    res = tabulate(k=2, r=1, max_stoich=2)
    assert res["total"] == 8
    assert res["by_deficiency"] == {0: 8}
    assert res["weakly_reversible_count"] == 0


def test_tabulate_k2_r2():
    res = tabulate(k=2, r=2, max_stoich=2)
    assert res["total"] == 99
    assert res["by_deficiency"] == {0: 82, 1: 17}
    assert res["weakly_reversible_count"] == 5


def test_tabulate_k3_r1():
    res = tabulate(k=3, r=1, max_stoich=2)
    assert res["total"] == 5
    assert res["by_deficiency"] == {0: 5}


@pytest.mark.slow
def test_tabulate_k2_r3():
    res = tabulate(k=2, r=3, max_stoich=2)
    assert res["total"] == 570
    assert res["by_deficiency"] == {0: 110, 1: 448, 2: 12}


@pytest.mark.slow
def test_tabulate_k3_r2():
    res = tabulate(k=3, r=2, max_stoich=2)
    assert res["total"] == 350
    assert res["by_deficiency"] == {0: 339, 1: 11}
