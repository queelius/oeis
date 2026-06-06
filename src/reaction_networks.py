"""
Chemical reaction networks: deficiency, weak reversibility, enumeration.

A *reaction network* on $k$ species is a finite set of directed reactions
$y \\to y'$ where each *complex* $y$ is a multiset over the species,
represented as a vector in $\\mathbb{Z}_{\\ge 0}^k$.  We work with networks
satisfying:

  - Reactant complex differs from product complex on each reaction.
  - Each complex has total stoichiometry (sum of entries) at most a bound
    (default 2: only mono- and bi-molecular reactions).
  - The reaction set is irreducible (no duplicate reactions).

Key invariants (Feinberg, "Foundations of Chemical Reaction Network Theory",
Springer 2019):

  - $n$ = number of distinct complexes appearing in the network.
  - $\\ell$ = number of *linkage classes* = weakly connected components of
    the directed graph whose vertices are complexes and whose edges are
    the reactions (forgetting direction).
  - $s$ = rank over $\\mathbb{Q}$ of the *stoichiometry matrix* whose
    columns are $y' - y$ for each reaction.
  - $\\delta = n - \\ell - s \\ge 0$ is the *deficiency*.
  - The network is *weakly reversible* if every reaction lies on a directed
    cycle in the reaction graph.

Feinberg's *Deficiency Zero Theorem*: if $\\delta = 0$ and the network is
weakly reversible, then for every choice of (positive) rate constants the
mass-action ODE admits a unique positive steady state, locally
asymptotically stable inside each stoichiometric compatibility class.

This module provides:

  - ``Network``: dataclass for a reaction network.
  - ``deficiency(net)``: integer $\\delta$.
  - ``linkage_classes(net)``: weakly connected components of complexes.
  - ``stoichiometry_rank(net)``: $s$ via sympy exact rational rank.
  - ``is_weakly_reversible(net)``: True iff every reaction is on a cycle.
  - ``is_deficiency_zero(net)``: True iff $\\delta = 0$.
  - ``canonical_form(net)``: relabel species lex-min, sort reactions.
  - ``enumerate_networks(k, r, max_stoich)``: enumerate all networks (up
    to canonical equivalence) with $k$ species, exactly $r$ reactions,
    and complexes of total stoichiometry $\\le$ max_stoich.

References:
  - Feinberg, M. (2019). Foundations of Chemical Reaction Network Theory.
    Applied Mathematical Sciences 202, Springer.
  - Feinberg, M. (1987). Chemical reaction network structure and the
    stability of complex isothermal reactors -- I. Chem. Eng. Sci. 42, 2229-2268.
  - Horn, F. and Jackson, R. (1972). General mass action kinetics.
    Arch. Rational Mech. Anal. 47, 81-116.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, permutations, product
from typing import Iterable

from sympy import Matrix


Complex = tuple[int, ...]  # one entry per species, count of that species
Reaction = tuple[Complex, Complex]  # (reactant, product)


# ---------------------------------------------------------------------------
# Core data class and invariants
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Network:
    """A chemical reaction network on a fixed number of species.

    Attributes:
        k: number of species.
        reactions: tuple of (reactant_complex, product_complex) pairs, each
            complex a tuple of length ``k`` with non-negative integer entries.
            All complexes implicitly defined by appearing in some reaction.
    """

    k: int
    reactions: tuple[Reaction, ...]

    def complexes(self) -> tuple[Complex, ...]:
        """Return distinct complexes appearing in any reaction, sorted lex."""
        seen: set[Complex] = set()
        for y, yp in self.reactions:
            seen.add(y)
            seen.add(yp)
        return tuple(sorted(seen))


def _complex_total(c: Complex) -> int:
    return sum(c)


def stoichiometry_rank(net: Network) -> int:
    """Rank over Q of the matrix whose columns are y' - y for each reaction.

    Uses sympy ``Matrix.rank()`` which is exact over the rationals.
    """
    if not net.reactions:
        return 0
    cols = []
    for y, yp in net.reactions:
        cols.append([yp[i] - y[i] for i in range(net.k)])
    # Matrix takes rows; we want columns to be reaction vectors, so transpose.
    M = Matrix(cols).T
    return int(M.rank())


def linkage_classes(net: Network) -> list[set[Complex]]:
    """Weakly connected components of the reaction graph.

    Vertices are the complexes; an undirected edge connects $y$ and $y'$
    whenever there is a reaction $y \\to y'$ or $y' \\to y$.
    """
    parent: dict[Complex, Complex] = {}

    def find(x: Complex) -> Complex:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: Complex, b: Complex) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for c in net.complexes():
        parent[c] = c
    for y, yp in net.reactions:
        union(y, yp)

    groups: dict[Complex, set[Complex]] = {}
    for c in net.complexes():
        r = find(c)
        groups.setdefault(r, set()).add(c)
    return list(groups.values())


def deficiency(net: Network) -> int:
    """Deficiency $\\delta = n - \\ell - s$ of the network.

    Always $\\ge 0$ (Feinberg).  We check this and raise if violated, since a
    negative value indicates a bug.
    """
    n = len(net.complexes())
    ell = len(linkage_classes(net))
    s = stoichiometry_rank(net)
    delta = n - ell - s
    if delta < 0:
        raise ValueError(
            f"Computed negative deficiency for {net!r}: "
            f"n={n}, ell={ell}, s={s}, delta={delta}. "
            "This contradicts Feinberg's theorem -- bug in the implementation."
        )
    return delta


def is_deficiency_zero(net: Network) -> bool:
    return deficiency(net) == 0


def is_weakly_reversible(net: Network) -> bool:
    """True iff every reaction lies on a directed cycle of the reaction graph.

    Equivalent: every weakly connected component is strongly connected.
    """
    if not net.reactions:
        return True
    # Build directed adjacency on complexes.
    adj: dict[Complex, set[Complex]] = {}
    for c in net.complexes():
        adj[c] = set()
    for y, yp in net.reactions:
        adj[y].add(yp)

    # Tarjan SCC.
    index_counter = [0]
    stack: list[Complex] = []
    on_stack: dict[Complex, bool] = {c: False for c in adj}
    indices: dict[Complex, int] = {}
    lowlinks: dict[Complex, int] = {}
    scc_id: dict[Complex, int] = {}
    next_scc = [0]

    def strongconnect(v: Complex) -> None:
        indices[v] = index_counter[0]
        lowlinks[v] = index_counter[0]
        index_counter[0] += 1
        stack.append(v)
        on_stack[v] = True
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                lowlinks[v] = min(lowlinks[v], lowlinks[w])
            elif on_stack[w]:
                lowlinks[v] = min(lowlinks[v], indices[w])
        if lowlinks[v] == indices[v]:
            sid = next_scc[0]
            next_scc[0] += 1
            while True:
                w = stack.pop()
                on_stack[w] = False
                scc_id[w] = sid
                if w == v:
                    break

    import sys

    sys.setrecursionlimit(10_000)
    for v in adj:
        if v not in indices:
            strongconnect(v)

    # Every reaction y -> yp must satisfy scc_id[y] == scc_id[yp].
    for y, yp in net.reactions:
        if scc_id[y] != scc_id[yp]:
            return False
    return True


# ---------------------------------------------------------------------------
# Canonical form
# ---------------------------------------------------------------------------


def _relabel(net: Network, perm: tuple[int, ...]) -> Network:
    """Apply a species permutation perm to all complexes.

    perm[i] = new index of old species i.
    """
    k = net.k
    new_reactions = []
    for y, yp in net.reactions:
        ny = [0] * k
        nyp = [0] * k
        for i in range(k):
            ny[perm[i]] = y[i]
            nyp[perm[i]] = yp[i]
        new_reactions.append((tuple(ny), tuple(nyp)))
    return Network(k=k, reactions=tuple(sorted(set(new_reactions))))


def _signature(net: Network) -> tuple:
    """Lex signature: sorted reactions as a tuple."""
    return tuple(sorted(net.reactions))


def canonical_form(net: Network) -> Network:
    """Return the lex-minimal network in the orbit under species permutation.

    With $k$ species the orbit has size at most $k!$, so for small $k$ this
    is cheap.  Two networks are isomorphic (via species relabeling) iff
    their canonical forms are equal.
    """
    best = None
    for perm in permutations(range(net.k)):
        cand = _relabel(net, perm)
        sig = _signature(cand)
        if best is None or sig < _signature(best):
            best = cand
    assert best is not None
    return best


# ---------------------------------------------------------------------------
# Enumeration of small networks
# ---------------------------------------------------------------------------


def enumerate_complexes(k: int, max_stoich: int) -> list[Complex]:
    """All complexes on $k$ species with total stoichiometry $\\le$ max_stoich.

    Enumerated lex-sorted; includes the zero complex (the "empty" complex).
    """
    out: list[Complex] = []
    # Complex c is a vector in Z_{>=0}^k with sum(c) in [0, max_stoich].
    def rec(prefix: list[int], remaining_slots: int, remaining_total: int) -> None:
        if remaining_slots == 0:
            out.append(tuple(prefix))
            return
        for v in range(remaining_total + 1):
            prefix.append(v)
            rec(prefix, remaining_slots - 1, remaining_total - v)
            prefix.pop()

    rec([], k, max_stoich)
    return sorted(out)


def enumerate_reactions(
    k: int, max_stoich: int, allow_zero_complex: bool = True
) -> list[Reaction]:
    """All directed reactions y -> yp with y != yp, both complexes within bound.

    If allow_zero_complex is False, exclude the zero complex (which
    corresponds to in/out flow from outside the system).
    """
    complexes = enumerate_complexes(k, max_stoich)
    if not allow_zero_complex:
        complexes = [c for c in complexes if any(c)]
    out: list[Reaction] = []
    for y in complexes:
        for yp in complexes:
            if y != yp:
                out.append((y, yp))
    return out


def enumerate_networks(
    k: int,
    r: int,
    max_stoich: int = 2,
    allow_zero_complex: bool = False,
    require_all_species: bool = True,
) -> Iterable[Network]:
    """Enumerate all networks on $k$ species with exactly $r$ reactions.

    Yields each isomorphism class exactly once (up to species permutation).
    Each complex has total stoichiometry $\\le$ max_stoich.

    Args:
        k: number of species.
        r: number of reactions.
        max_stoich: max sum of complex entries (default 2 = bimolecular).
        allow_zero_complex: include the all-zeros complex (in/out flows)?
        require_all_species: only keep networks in which every species
            appears in some complex.  If False, networks that are really
            $j$-species networks for $j < k$ are also included.

    This is an unfiltered Cartesian-product-then-deduplicate approach, only
    tractable for small $(k, r, \\text{max\\_stoich})$.  For $k=2, r\\le 4,
    \\text{max\\_stoich}=2$ the search is very small; for larger parameters
    one should pre-filter by isomorphism during enumeration.
    """
    reactions = enumerate_reactions(k, max_stoich, allow_zero_complex=allow_zero_complex)
    seen: set[tuple] = set()
    for combo in combinations(reactions, r):
        net = Network(k=k, reactions=tuple(sorted(combo)))
        if require_all_species:
            used = set()
            for y, yp in net.reactions:
                for i in range(k):
                    if y[i] or yp[i]:
                        used.add(i)
            if len(used) != k:
                continue
        canon = canonical_form(net)
        sig = _signature(canon)
        if sig in seen:
            continue
        seen.add(sig)
        yield canon


# ---------------------------------------------------------------------------
# Tabulation helpers
# ---------------------------------------------------------------------------


def tabulate(
    k: int,
    r: int,
    max_stoich: int = 2,
    allow_zero_complex: bool = False,
    require_all_species: bool = True,
) -> dict:
    """Tabulate (deficiency, weakly_reversible) distribution for given (k, r).

    Returns a dict with:
        total: total number of inequivalent networks.
        by_deficiency: mapping delta -> count.
        by_wr: mapping (delta, is_wr) -> count.
        max_deficiency: largest delta seen.
        weakly_reversible_count: total number with weak reversibility.
    """
    nets = list(
        enumerate_networks(
            k=k,
            r=r,
            max_stoich=max_stoich,
            allow_zero_complex=allow_zero_complex,
            require_all_species=require_all_species,
        )
    )
    by_def: dict[int, int] = {}
    by_wr: dict[tuple[int, bool], int] = {}
    wr_count = 0
    for net in nets:
        d = deficiency(net)
        by_def[d] = by_def.get(d, 0) + 1
        wr = is_weakly_reversible(net)
        if wr:
            wr_count += 1
        by_wr[(d, wr)] = by_wr.get((d, wr), 0) + 1
    return {
        "total": len(nets),
        "by_deficiency": dict(sorted(by_def.items())),
        "by_wr": dict(sorted(by_wr.items())),
        "max_deficiency": max(by_def) if by_def else 0,
        "weakly_reversible_count": wr_count,
    }


# ---------------------------------------------------------------------------
# Pretty printing for sanity checks
# ---------------------------------------------------------------------------


def format_complex(c: Complex, names: tuple[str, ...] | None = None) -> str:
    if names is None:
        names = tuple(chr(ord("A") + i) for i in range(len(c)))
    parts: list[str] = []
    for i, v in enumerate(c):
        if v == 0:
            continue
        if v == 1:
            parts.append(names[i])
        else:
            parts.append(f"{v}{names[i]}")
    return "0" if not parts else "+".join(parts)


def format_network(net: Network, names: tuple[str, ...] | None = None) -> str:
    if names is None:
        names = tuple(chr(ord("A") + i) for i in range(net.k))
    lines = []
    for y, yp in net.reactions:
        lines.append(f"  {format_complex(y, names)} -> {format_complex(yp, names)}")
    return "\n".join(lines)
