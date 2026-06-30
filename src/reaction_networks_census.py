"""Complete census of small chemical reaction networks by Feinberg deficiency.

This module *extends* :mod:`src.reaction_networks`.  Where that module
tabulates networks with **exactly** ``r`` reactions, here we build the
**complete census**: every reaction network (up to species-permutation
``S_k`` canonicalization) whose complexes are mono- or bi-molecular
(total stoichiometry ``<= max_stoich``, default 2), with no zero complex,
and with every species used.  The census is bounded -- the ground set of
directed reactions is finite -- so the deficiency distribution over *all*
networks (any number of reactions) is a well-defined finite object.

Census objects
--------------

For fixed ``(k, max_stoich)`` let ``R`` be the number of directed
reactions ``y -> y'`` with ``y != y'``, both complexes non-zero with
total stoichiometry ``<= max_stoich``.  A *network* is a non-empty subset
of these ``R`` reactions in which every species appears.  Two networks
are identified iff a species permutation in ``S_k`` carries one to the
other.  We report, over all such networks:

  - the total count,
  - the count by deficiency ``delta = n - ell - s``,
  - the count of weakly reversible networks,
  - the count of weakly reversible *and* deficiency-zero networks
    (the hypothesis class of Feinberg's Deficiency Zero Theorem),
  - the same distributions refined by the number of reactions ``r``.

Feasibility
-----------

The ground set has size ``R``; a complete census enumerates the ``2^R``
subsets.

  ============  ==============  ====================
  ``k``          ``R``           complete census?
  ============  ==============  ====================
  1              2               yes (``2^2``)
  2              20              yes (``2^20 ~ 1e6``)
  3              72              no (``2^72``); bounded by ``r``
  4              182             no; bounded by small ``r``
  ============  ==============  ====================

For ``k <= 2`` we perform the *complete* census (every subset).  For
``k >= 3`` the complete census is infeasible, so we enumerate by reaction
count ``r`` up to a bound (``C(R, r)`` subsets), which is what the
original module already does; this module adds a faster exact-rank path
and a Burnside cross-check of the *labelled* totals.

Performance
-----------

The two hot paths from :mod:`src.reaction_networks` are replaced:

  - ``stoichiometry_rank`` (sympy) -> :func:`integer_rank`, a fraction-free
    (Bareiss) integer elimination, validated bit-for-bit against sympy.
  - canonicalization under ``S_k`` is done inline on the frozen reaction
    set; for ``k <= 2`` only ``k!`` permutations are tried.

The ``k = 2`` complete census (``2^20`` subsets) is parallelised across
processes via :func:`census` with ``workers > 1``.

References
----------
  - Feinberg, M. (2019). Foundations of Chemical Reaction Network Theory.
    Applied Mathematical Sciences 202, Springer.  (Deficiency, the
    Deficiency Zero Theorem.)
  - Horn, F. and Jackson, R. (1972). General mass action kinetics.
    Arch. Rational Mech. Anal. 47, 81-116.
"""

from __future__ import annotations

import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import combinations, permutations
from math import factorial
from typing import Iterable

from src.reaction_networks import (
    Complex,
    Network,
    Reaction,
    enumerate_reactions,
)

# ---------------------------------------------------------------------------
# Fast exact stoichiometric rank (fraction-free, integer arithmetic)
# ---------------------------------------------------------------------------


def integer_rank(rows: list[list[int]]) -> int:
    """Rank over Q of an integer matrix via fraction-free Bareiss elimination.

    Operates on Python ints (arbitrary precision), so the result is exact.
    ``rows`` is a list of equal-length integer rows.  Returns the rank.

    The Bareiss algorithm keeps all intermediate entries integral, avoiding
    both floating point and rational arithmetic.  Validated bit-for-bit
    against ``sympy.Matrix.rank`` in the test suite.
    """
    if not rows:
        return 0
    M = [list(map(int, r)) for r in rows]
    nrows = len(M)
    ncols = len(M[0])
    rank = 0
    prev_pivot = 1
    pivot_row = 0
    for col in range(ncols):
        # locate a non-zero pivot at or below pivot_row in this column
        sel = -1
        for i in range(pivot_row, nrows):
            if M[i][col] != 0:
                sel = i
                break
        if sel == -1:
            continue
        if sel != pivot_row:
            M[pivot_row], M[sel] = M[sel], M[pivot_row]
        pivot_val = M[pivot_row][col]
        for i in range(nrows):
            if i == pivot_row:
                continue
            # fraction-free update: M[i][j] = (pivot*M[i][j] - M[i][col]*M[p][j]) / prev_pivot
            factor = M[i][col]
            if factor == 0:
                # still must update for the Bareiss divisor bookkeeping below;
                # but with factor 0 the row is unchanged up to the /prev_pivot,
                # which for an all-integer matrix stays integral.
                row_i = M[i]
                row_p = M[pivot_row]
                for j in range(col, ncols):
                    row_i[j] = (pivot_val * row_i[j]) // prev_pivot
                continue
            row_i = M[i]
            row_p = M[pivot_row]
            for j in range(col, ncols):
                row_i[j] = (pivot_val * row_i[j] - factor * row_p[j]) // prev_pivot
        prev_pivot = pivot_val
        pivot_row += 1
        rank += 1
        if pivot_row == nrows:
            break
    return rank


# ---------------------------------------------------------------------------
# Fast invariants on a raw reaction tuple (no Network allocation overhead)
# ---------------------------------------------------------------------------


def _complexes_of(reactions: tuple[Reaction, ...]) -> list[Complex]:
    seen: dict[Complex, None] = {}
    for y, yp in reactions:
        seen.setdefault(y, None)
        seen.setdefault(yp, None)
    return list(seen)


def _num_linkage_classes(reactions: tuple[Reaction, ...], complexes: list[Complex]) -> int:
    """Weakly connected components of the reaction graph (union-find)."""
    parent: dict[Complex, Complex] = {c: c for c in complexes}

    def find(x: Complex) -> Complex:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for y, yp in reactions:
        ra, rb = find(y), find(yp)
        if ra != rb:
            parent[ra] = rb
    roots = {find(c) for c in complexes}
    return len(roots)


def _stoich_rank(reactions: tuple[Reaction, ...], k: int) -> int:
    rows = [[yp[i] - y[i] for i in range(k)] for (y, yp) in reactions]
    return integer_rank(rows)


def deficiency_fast(reactions: tuple[Reaction, ...], k: int) -> int:
    """Deficiency ``n - ell - s`` computed without allocating a Network."""
    complexes = _complexes_of(reactions)
    n = len(complexes)
    ell = _num_linkage_classes(reactions, complexes)
    s = _stoich_rank(reactions, k)
    delta = n - ell - s
    if delta < 0:
        raise ValueError(
            f"negative deficiency n={n} ell={ell} s={s} for {reactions!r}"
        )
    return delta


def is_weakly_reversible_fast(reactions: tuple[Reaction, ...]) -> bool:
    """True iff every reaction lies on a directed cycle (every WCC is an SCC).

    Iterative Tarjan to avoid Python recursion limits on larger complex sets.
    """
    if not reactions:
        return True
    adj: dict[Complex, list[Complex]] = {}
    for y, yp in reactions:
        adj.setdefault(y, [])
        adj.setdefault(yp, [])
    for y, yp in reactions:
        adj[y].append(yp)

    index = 0
    indices: dict[Complex, int] = {}
    lowlink: dict[Complex, int] = {}
    on_stack: dict[Complex, bool] = {}
    stack: list[Complex] = []
    scc_id: dict[Complex, int] = {}
    next_scc = 0

    for start in adj:
        if start in indices:
            continue
        # iterative DFS: work items are (node, neighbor_iterator_position)
        work: list[tuple[Complex, int]] = [(start, 0)]
        while work:
            v, pi = work[-1]
            if pi == 0:
                indices[v] = index
                lowlink[v] = index
                index += 1
                stack.append(v)
                on_stack[v] = True
            recursed = False
            neighbors = adj[v]
            while pi < len(neighbors):
                w = neighbors[pi]
                pi += 1
                if w not in indices:
                    work[-1] = (v, pi)
                    work.append((w, 0))
                    recursed = True
                    break
                elif on_stack.get(w, False):
                    if indices[w] < lowlink[v]:
                        lowlink[v] = indices[w]
            if recursed:
                continue
            work[-1] = (v, pi)
            # done exploring v's neighbors
            if lowlink[v] == indices[v]:
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    scc_id[w] = next_scc
                    if w == v:
                        break
                next_scc += 1
            work.pop()
            if work:
                parent_v = work[-1][0]
                if lowlink[v] < lowlink[parent_v]:
                    lowlink[parent_v] = lowlink[v]

    for y, yp in reactions:
        if scc_id[y] != scc_id[yp]:
            return False
    return True


# ---------------------------------------------------------------------------
# Canonicalization under S_k (inline, fast)
# ---------------------------------------------------------------------------


def _perms(k: int) -> list[tuple[int, ...]]:
    return list(permutations(range(k)))


def canonical_signature(
    reactions: tuple[Reaction, ...], k: int, perms: list[tuple[int, ...]]
) -> tuple[Reaction, ...]:
    """Lex-min reaction tuple over all species permutations in ``perms``.

    Returns the canonical (sorted) reaction tuple of the orbit, suitable as
    a hashable key identifying the isomorphism class.
    """
    best: tuple[Reaction, ...] | None = None
    for perm in perms:
        relabeled = []
        for y, yp in reactions:
            ny = [0] * k
            nyp = [0] * k
            for i in range(k):
                ny[perm[i]] = y[i]
                nyp[perm[i]] = yp[i]
            relabeled.append((tuple(ny), tuple(nyp)))
        cand = tuple(sorted(relabeled))
        if best is None or cand < best:
            best = cand
    assert best is not None
    return best


def _all_species_used(reactions: tuple[Reaction, ...], k: int) -> bool:
    used = 0
    for y, yp in reactions:
        for i in range(k):
            if y[i] or yp[i]:
                used |= 1 << i
    return used == (1 << k) - 1


# ---------------------------------------------------------------------------
# Census result container
# ---------------------------------------------------------------------------


@dataclass
class CensusResult:
    """Aggregated deficiency distribution of a network census."""

    k: int
    max_stoich: int
    complete: bool  # True if the full 2^R census; False if bounded by r
    r_max: int | None  # reaction-count bound when not complete
    total: int
    by_deficiency: dict[int, int]
    weakly_reversible: int
    wr_and_def0: int
    # refinement by reaction count r: r -> {"total", "by_deficiency", "wr", "wr_def0"}
    by_r: dict[int, dict]

    def deficiency_sequence(self) -> list[int]:
        """Counts ``[c_0, c_1, ...]`` indexed by deficiency value."""
        if not self.by_deficiency:
            return []
        m = max(self.by_deficiency)
        return [self.by_deficiency.get(d, 0) for d in range(m + 1)]

    def r_total_sequence(self) -> list[int]:
        """Total network count indexed by reaction count r=1,2,..."""
        if not self.by_r:
            return []
        m = max(self.by_r)
        return [self.by_r.get(r, {}).get("total", 0) for r in range(1, m + 1)]


# ---------------------------------------------------------------------------
# Subset enumeration
# ---------------------------------------------------------------------------


def _process_chunk(
    args: tuple[list[Reaction], int, list[tuple[int, ...]], int, int, int]
) -> dict:
    """Worker: process bitmasks in ``[lo, hi)`` over the reaction ground set.

    Returns partial canonical-form sets and tallies.  Because different
    workers may independently discover the same canonical form, we return the
    *set of canonical signatures* (deduplicated in the parent) together with,
    for each signature, its (delta, wr, r) -- these are functions of the
    canonical form alone, so duplicates agree and can be merged by union.
    """
    reactions, k, perms, R, lo, hi = args
    # local map signature -> (delta, wr, r)
    local: dict[tuple, tuple[int, bool, int]] = {}
    full_mask = (1 << k) - 1
    for mask in range(lo, hi):
        if mask == 0:
            continue
        # all species used?  cheap pre-filter on the mask's reactions
        combo = tuple(reactions[i] for i in range(R) if (mask >> i) & 1)
        used = 0
        for y, yp in combo:
            for i in range(k):
                if y[i] or yp[i]:
                    used |= 1 << i
            if used == full_mask:
                break
        if used != full_mask:
            continue
        sig = canonical_signature(combo, k, perms)
        if sig in local:
            continue
        delta = deficiency_fast(sig, k)
        wr = is_weakly_reversible_fast(sig)
        local[sig] = (delta, wr, len(sig))
    return local


def census(
    k: int,
    max_stoich: int = 2,
    r_max: int | None = None,
    allow_zero_complex: bool = False,
    workers: int | None = None,
) -> CensusResult:
    """Complete (or r-bounded) census of networks by deficiency.

    Args:
        k: number of species.
        max_stoich: max total stoichiometry per complex (default 2).
        r_max: if given, only enumerate networks with at most this many
            reactions (bounds the work to ``sum_{r<=r_max} C(R, r)``).  If
            ``None``, enumerate the complete ``2^R`` census; only feasible
            for small ``k`` (k <= 2 with max_stoich=2).
        allow_zero_complex: include the zero complex (in/out flows)?
        workers: process pool size for the complete census.  Defaults to all
            cores when the full census is large, else single-threaded.

    Returns:
        A :class:`CensusResult`.
    """
    reactions = enumerate_reactions(k, max_stoich, allow_zero_complex=allow_zero_complex)
    R = len(reactions)
    perms = _perms(k)

    # decide enumeration strategy
    complete = r_max is None
    if complete:
        return _census_complete(k, max_stoich, reactions, R, perms, workers)
    if workers is not None and workers > 1:
        return _census_by_r_parallel(
            k, max_stoich, reactions, R, perms, r_max, workers
        )
    return _census_by_r(k, max_stoich, reactions, R, perms, r_max)


def _aggregate(local: dict[tuple, tuple[int, bool, int]], k: int, max_stoich: int,
               complete: bool, r_max: int | None) -> CensusResult:
    by_def: Counter = Counter()
    wr_total = 0
    wr_def0 = 0
    by_r: dict[int, dict] = {}
    for sig, (delta, wr, r) in local.items():
        by_def[delta] += 1
        if wr:
            wr_total += 1
            if delta == 0:
                wr_def0 += 1
        slot = by_r.setdefault(
            r, {"total": 0, "by_deficiency": Counter(), "wr": 0, "wr_def0": 0}
        )
        slot["total"] += 1
        slot["by_deficiency"][delta] += 1
        if wr:
            slot["wr"] += 1
            if delta == 0:
                slot["wr_def0"] += 1
    # freeze counters into dicts
    for r in by_r:
        by_r[r]["by_deficiency"] = dict(sorted(by_r[r]["by_deficiency"].items()))
    return CensusResult(
        k=k,
        max_stoich=max_stoich,
        complete=complete,
        r_max=r_max,
        total=len(local),
        by_deficiency=dict(sorted(by_def.items())),
        weakly_reversible=wr_total,
        wr_and_def0=wr_def0,
        by_r=dict(sorted(by_r.items())),
    )


def _census_complete(
    k: int,
    max_stoich: int,
    reactions: list[Reaction],
    R: int,
    perms: list[tuple[int, ...]],
    workers: int | None,
) -> CensusResult:
    n_subsets = 1 << R
    if workers is None:
        workers = (os.cpu_count() or 1) if n_subsets >= (1 << 16) else 1

    if workers <= 1 or n_subsets < 1 << 14:
        local = _process_chunk((reactions, k, perms, R, 0, n_subsets))
        return _aggregate(local, k, max_stoich, True, None)

    # split [0, 2^R) into chunks across workers
    n_chunks = workers * 4
    step = (n_subsets + n_chunks - 1) // n_chunks
    tasks = []
    lo = 0
    while lo < n_subsets:
        hi = min(lo + step, n_subsets)
        tasks.append((reactions, k, perms, R, lo, hi))
        lo = hi

    merged: dict[tuple, tuple[int, bool, int]] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for partial in ex.map(_process_chunk, tasks):
            merged.update(partial)
    return _aggregate(merged, k, max_stoich, True, None)


def _census_by_r(
    k: int,
    max_stoich: int,
    reactions: list[Reaction],
    R: int,
    perms: list[tuple[int, ...]],
    r_max: int,
) -> CensusResult:
    local: dict[tuple, tuple[int, bool, int]] = {}
    full_mask = (1 << k) - 1
    for r in range(1, r_max + 1):
        for combo in combinations(reactions, r):
            used = 0
            for y, yp in combo:
                for i in range(k):
                    if y[i] or yp[i]:
                        used |= 1 << i
            if used != full_mask:
                continue
            sig = canonical_signature(combo, k, perms)
            if sig in local:
                continue
            delta = deficiency_fast(sig, k)
            wr = is_weakly_reversible_fast(sig)
            local[sig] = (delta, wr, len(sig))
    return _aggregate(local, k, max_stoich, False, r_max)


def _process_combo_list(
    args: tuple[tuple[Reaction, ...], int, list[tuple[int, ...]], int, list]
) -> dict:
    """Worker for the parallel by-r census over an explicit list of combos.

    The parent slices ``combinations(reactions, r)`` *sequentially* (so the
    underlying iterator is advanced exactly once, with no per-worker re-skip)
    and ships each materialized chunk here.  Returns the local
    signature -> (delta, wr, r) map; the parent merges by union (each value is
    a function of the canonical form, so collisions agree).
    """
    reactions, k, perms, r, combos = args
    local: dict[tuple, tuple[int, bool, int]] = {}
    full_mask = (1 << k) - 1
    for combo in combos:
        used = 0
        for y, yp in combo:
            for i in range(k):
                if y[i] or yp[i]:
                    used |= 1 << i
            if used == full_mask:
                break
        if used != full_mask:
            continue
        sig = canonical_signature(combo, k, perms)
        if sig in local:
            continue
        delta = deficiency_fast(sig, k)
        wr = is_weakly_reversible_fast(sig)
        local[sig] = (delta, wr, len(sig))
    return local


def _census_by_r_parallel(
    k: int,
    max_stoich: int,
    reactions: list[Reaction],
    R: int,
    perms: list[tuple[int, ...]],
    r_max: int,
    workers: int,
    chunk_size: int = 200_000,
) -> CensusResult:
    """Parallel by-r census.

    The combinations iterator is advanced exactly once in the parent and cut
    into contiguous chunks of ``chunk_size`` combos via sequential
    :func:`itertools.islice` (each ``islice`` resumes where the previous left
    off -- no quadratic re-skip).  Chunks are dispatched to the pool with
    bounded look-ahead so memory stays flat regardless of C(R, r).
    """
    from itertools import islice

    merged: dict[tuple, tuple[int, bool, int]] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in range(1, r_max + 1):
            it = combinations(reactions, r)
            # bounded-lookahead dispatch: keep at most 2*workers chunks in flight
            from collections import deque

            inflight: deque = deque()
            max_inflight = 2 * workers

            def submit_next() -> bool:
                chunk = list(islice(it, chunk_size))
                if not chunk:
                    return False
                inflight.append(
                    ex.submit(
                        _process_combo_list, (reactions, k, perms, r, chunk)
                    )
                )
                return True

            for _ in range(max_inflight):
                if not submit_next():
                    break
            while inflight:
                fut = inflight.popleft()
                merged.update(fut.result())
                submit_next()
    return _aggregate(merged, k, max_stoich, False, r_max)


# ---------------------------------------------------------------------------
# Burnside cross-check on labelled totals
# ---------------------------------------------------------------------------


def labelled_total(k: int, max_stoich: int, r: int,
                   allow_zero_complex: bool = False,
                   require_all_species: bool = True) -> int:
    """Number of *labelled* (not up to symmetry) networks with exactly r reactions.

    A labelled network is any r-subset of the R directed reactions in which
    every species is used (when ``require_all_species``).  This is an
    independent count used to cross-check the orbit count via Burnside.
    """
    reactions = enumerate_reactions(k, max_stoich, allow_zero_complex=allow_zero_complex)
    R = len(reactions)
    if not require_all_species:
        return _comb(R, r)
    full_mask = (1 << k) - 1
    count = 0
    for combo in combinations(reactions, r):
        used = 0
        for y, yp in combo:
            for i in range(k):
                if y[i] or yp[i]:
                    used |= 1 << i
        if used == full_mask:
            count += 1
    return count


def _comb(n: int, r: int) -> int:
    from math import comb

    return comb(n, r)


def orbit_count_burnside(k: int, max_stoich: int, r: int,
                         allow_zero_complex: bool = False) -> int:
    """Count r-reaction networks up to S_k via Burnside's lemma.

    The number of orbits equals the average over g in S_k of the number of
    g-fixed labelled r-subsets (with all species used).  Computed directly
    by, for each permutation g, counting r-subsets fixed by g.  This is an
    independent route to the orbit count produced by :func:`census` and the
    original ``tabulate``; it requires *all species used*, matching the
    census convention.

    Note: because "all species used" is itself an S_k-invariant predicate,
    Burnside applies to the restricted ground set of all-species-used subsets.
    """
    reactions = enumerate_reactions(k, max_stoich, allow_zero_complex=allow_zero_complex)
    R = len(reactions)
    full_mask = (1 << k) - 1
    perms = _perms(k)

    # action of each permutation on the reaction index set
    def relabel_reaction(rx: Reaction, perm: tuple[int, ...]) -> Reaction:
        y, yp = rx
        ny = [0] * k
        nyp = [0] * k
        for i in range(k):
            ny[perm[i]] = y[i]
            nyp[perm[i]] = yp[i]
        return (tuple(ny), tuple(nyp))

    index_of = {rx: i for i, rx in enumerate(reactions)}
    # species mask of each reaction (which species it touches)
    rmask = [0] * R
    for i, (y, yp) in enumerate(reactions):
        m = 0
        for j in range(k):
            if y[j] or yp[j]:
                m |= 1 << j
        rmask[i] = m

    total_fixed = 0
    for perm in perms:
        img = [index_of[relabel_reaction(reactions[i], perm)] for i in range(R)]
        # cycle decomposition of the action on reaction indices
        seen = [False] * R
        cycles: list[tuple[int, int]] = []  # (cycle_length, union_species_mask)
        for i in range(R):
            if seen[i]:
                continue
            length = 0
            cmask = 0
            j = i
            while not seen[j]:
                seen[j] = True
                length += 1
                cmask |= rmask[j]
                j = img[j]
            cycles.append((length, cmask))
        # A g-fixed subset is a union of whole cycles.  We need those of total
        # size r whose union species-mask is full_mask.  DP over cycles tracking
        # (accumulated species mask) -> {size -> count}.  Only 2^k masks and
        # sizes 0..r, so this is polynomial, not 2^{#cycles}.
        dp: list[list[int]] = [[0] * (r + 1) for _ in range(1 << k)]
        dp[0][0] = 1
        for clen, cmask in cycles:
            if clen > r:
                # taking this cycle can never reach exactly r if clen>r alone,
                # but it could still overshoot; skip-only contribution is the
                # identity, so just continue (cannot be part of a size-r union
                # unless combined to <= r, impossible since clen>r).
                continue
            new_dp = [row[:] for row in dp]  # skip case carried over
            for mask in range(1 << k):
                row = dp[mask]
                nmask = mask | cmask
                tgt = new_dp[nmask]
                for size in range(r - clen + 1):
                    v = row[size]
                    if v:
                        tgt[size + clen] += v
            dp = new_dp
        total_fixed += dp[full_mask][r]
    assert total_fixed % len(perms) == 0
    return total_fixed // len(perms)
