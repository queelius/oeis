#!/usr/bin/env python3
"""
Rainbow (anti-Ramsey) Schur numbers for finite abelian groups.

For a finite abelian group G and the Schur equation x1 + x2 = x3, the
rainbow number rb(G) is the smallest number of colors r such that EVERY
exact (surjective) r-coloring of the ground set admits a *rainbow* Schur
triple: a set {x1, x2, x3} of three group elements with x1 + x2 = x3 that
receives three pairwise-distinct colors.

Equivalently,

    rb(G) = 1 + (max number of colors in a rainbow-free exact coloring).

Convention (locked against the cyclic and grid literature):
  * coloring c : S -> [r], "exact" means surjective (every color used);
  * a "solution" is a SET {x1, x2, x3} with x1 + x2 = x3;
  * a solution with a repeated element is *degenerate* and CANNOT be
    rainbow (e.g. 2x = z, i.e. x1 = x2, is never rainbow, since
    c(x1) = c(x2)); equivalently a rainbow triple needs x1, x2, x3
    pairwise distinct as elements AND pairwise distinct in color;
  * if S has no (non-degenerate) solution at all, rb(G) = |S| + 1.

Ground-set convention is a parameter:
  * include_zero=True  : color ALL elements of G  (the convention used by
    Bevilacqua-King-Kritschgau-Tait-Tebon-Young, arXiv:1809.04576, for
    cyclic Z_n: c : Z_n -> [r]);
  * include_zero=False : color only the nonzero elements G minus {0}.
The validation gate below establishes which convention reproduces the
published cyclic values rb(Z_p,1)=4 (p>=5 prime), rb(Z_2)=rb(Z_3)=3.

References
----------
Bevilacqua, King, Kritschgau, Tait, Tebon, Young.
  "Rainbow numbers for x1 + x2 = k x3 in Z_n." arXiv:1809.04576.
Fallon, Manhart, Miller, Rehm, Warnberg, Zinnel.
  "Rainbow numbers of [m] x [n] for x1 + x2 = x3." arXiv:2301.10349.

This module is pure standard library (plus an optional python-sat backend
for the scalable solver). The group machinery mirrors src/schur_groups.py.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

Elem = Tuple[int, ...]


# --------------------------------------------------------------------------
# Group construction: G = Z_{n1} x Z_{n2} x ... as tuples, componentwise add
# --------------------------------------------------------------------------

def group_elements(orders: Sequence[int]) -> List[Elem]:
    """All elements of Z_{n1} x ... x Z_{nk} as tuples, identity first-ish."""
    return list(itertools.product(*[range(n) for n in orders]))


def group_add(a: Elem, b: Elem, orders: Sequence[int]) -> Elem:
    return tuple((x + y) % n for x, y, n in zip(a, b, orders))


def group_zero(orders: Sequence[int]) -> Elem:
    return tuple(0 for _ in orders)


def group_order(orders: Sequence[int]) -> int:
    o = 1
    for n in orders:
        o *= n
    return o


# --------------------------------------------------------------------------
# Schur triples
# --------------------------------------------------------------------------

def schur_triples(orders: Sequence[int], include_zero: bool = True
                  ) -> Tuple[List[Elem], List[Tuple[int, int, int]]]:
    """
    Ground set S and the list of non-degenerate Schur triples on it.

    Returns (elements, triples) where `elements` is the ground set (a list
    of group elements, indexed 0..|S|-1) and `triples` is a list of index
    triples (i, j, k) with i < j (unordered in x1, x2) such that
    elements[i] + elements[j] == elements[k] and i, j, k are pairwise
    distinct (degenerate triples are excluded -- they can never be rainbow).

    include_zero controls whether 0 is part of the colored ground set.
    """
    elems = group_elements(orders)
    if not include_zero:
        z = group_zero(orders)
        elems = [e for e in elems if e != z]
    index = {e: i for i, e in enumerate(elems)}
    n = len(elems)
    triples: List[Tuple[int, int, int]] = []
    for i in range(n):
        for j in range(i, n):
            s = group_add(elems[i], elems[j], orders)
            k = index.get(s)
            if k is None:
                continue
            # degenerate (repeated element) cannot be rainbow -> skip
            if i == j or k == i or k == j:
                continue
            triples.append((i, j, k))
    return elems, triples


# --------------------------------------------------------------------------
# Rainbow-freeness of a partition (a coloring given as a label per element)
# --------------------------------------------------------------------------

def is_rainbow_free(labels: Sequence[int],
                    triples: Sequence[Tuple[int, int, int]]) -> bool:
    """A labelling (color per ground element) is rainbow-free iff no Schur
    triple gets three distinct labels."""
    for i, j, k in triples:
        a, b, c = labels[i], labels[j], labels[k]
        if a != b and a != c and b != c:
            return False
    return True


# --------------------------------------------------------------------------
# Brute-force ground truth (set-partition enumeration), for small |S|
# --------------------------------------------------------------------------

def _set_partitions_by_blocks(n: int):
    """Yield (num_blocks, labels) for every set partition of {0,..,n-1}.

    Uses restricted-growth strings so labels are canonical (first element
    of each new block gets the next integer). This enumerates each
    *unordered* partition exactly once -- the right object, since rainbow-
    freeness and the color count depend only on the partition, not on the
    color names.
    """
    if n == 0:
        yield 0, []
        return
    a = [0] * n          # restricted growth string
    b = [0] * n          # b[i] = max(a[0..i-1]) i.e. running max+? see below
    # Knuth Algorithm H for restricted growth strings.
    # a[0]=0 always. b[i] = 1 + max(a[0..i]).
    # We implement the standard generator.
    a = [0] * n
    m = [0] * n  # m[i] = max label allowed at position i = 1 + max(a[:i])
    # iterative generation
    def gen(i, max_label):
        if i == n:
            yield max_label + 1, list(a)
            return
        for v in range(max_label + 1):
            a[i] = v
            yield from gen(i + 1, max_label)
        a[i] = max_label + 1
        yield from gen(i + 1, max_label + 1)

    # position 0 fixed to 0
    a[0] = 0
    yield from gen(1, 0)


def max_rainbow_free_colors_bruteforce(orders: Sequence[int],
                                       include_zero: bool = True) -> int:
    """Maximum number of colors in a rainbow-free exact coloring, by
    exhaustive enumeration of set partitions. Ground truth; |S| small."""
    elems, triples = schur_triples(orders, include_zero)
    n = len(elems)
    if n == 0:
        return 0
    best = 1
    for nb, labels in _set_partitions_by_blocks(n):
        if nb <= best:
            continue
        if is_rainbow_free(labels, triples):
            best = nb
    return best


def rb_bruteforce(orders: Sequence[int], include_zero: bool = True) -> int:
    """rb(G) by brute force. rb = (max rainbow-free colors) + 1, with the
    no-solution convention rb = |S| + 1."""
    elems, triples = schur_triples(orders, include_zero)
    n = len(elems)
    if n == 0:
        return 1  # empty ground set: vacuous, |S|+1 = 1
    if not triples:
        return n + 1  # no solutions at all -> convention
    return max_rainbow_free_colors_bruteforce(orders, include_zero) + 1


# --------------------------------------------------------------------------
# Scalable exact solver: "does a rainbow-free surjective r-coloring exist?"
# via SAT, searching r downward from |S| to find the maximum.
# --------------------------------------------------------------------------

def _exists_rainbow_free_with_r_colors_sat(
        n: int,
        triples: Sequence[Tuple[int, int, int]],
        r: int) -> bool:
    """SAT feasibility: exists a surjective r-coloring of an n-element set,
    rainbow-free w.r.t. `triples`.

    Encoding: boolean x[v,c] = "element v has color c", v in 0..n-1,
    c in 0..r-1.
      * each element gets >=1 color (ALO) and <=1 color (AMO);
      * each color used >=1 time (surjectivity);
      * symmetry break: element v may only use colors 0..v (so the first
        occurrence of color c is at position >= c); this kills the r!
        color-permutation symmetry and the colors-0..min(v,r-1) ordering;
      * rainbow-free: for every triple (i,j,k) and every choice of three
        DISTINCT colors (a,b,c) we forbid x[i,a]&x[j,b]&x[k,c]. Because the
        roles x1,x2 are symmetric (i<j) and x3 is k, we forbid all ordered
        distinct (a,b,c) assignments to (i,j,k).
    """
    from pysat.formula import CNF
    from pysat.solvers import Solver

    if r > n or r < 1:
        return False

    def var(v, c):  # 1-based DIMACS variable id
        return v * r + c + 1

    cnf = CNF()

    # Symmetry break for color labels: element v can use colors 0..min(v, r-1).
    # (first element fixes color 0; element v introduces at most one new color.)
    allowed = [list(range(min(v, r - 1) + 1)) for v in range(n)]

    # ALO + AMO per element (only over allowed colors)
    for v in range(n):
        cols = allowed[v]
        if not cols:
            return False
        cnf.append([var(v, c) for c in cols])         # ALO
        for a in range(len(cols)):                    # AMO pairwise
            for b in range(a + 1, len(cols)):
                cnf.append([-var(v, cols[a]), -var(v, cols[b])])
        # colors outside allowed are forced false (omit vars entirely by
        # never asserting them; but they could float -> pin to false)
        for c in range(r):
            if c not in cols:
                cnf.append([-var(v, c)])

    # Surjectivity: every color used at least once (only by elements allowed
    # to take it, i.e. v >= c).
    for c in range(r):
        clause = [var(v, c) for v in range(n) if c in allowed[v]]
        if not clause:
            return False
        cnf.append(clause)

    # Rainbow-free constraints.
    for (i, j, k) in triples:
        ci, cj, ck = allowed[i], allowed[j], allowed[k]
        for a in ci:
            for b in cj:
                if b == a:
                    continue
                for c in ck:
                    if c == a or c == b:
                        continue
                    cnf.append([-var(i, a), -var(j, b), -var(k, c)])

    with Solver(name="g3", bootstrap_with=cnf.clauses) as s:
        return s.solve()


def max_rainbow_free_colors_sat(orders: Sequence[int],
                                include_zero: bool = True,
                                upper: Optional[int] = None) -> int:
    """Maximum number of colors in a rainbow-free exact coloring, via SAT.
    Searches r downward from an upper bound until feasible."""
    elems, triples = schur_triples(orders, include_zero)
    n = len(elems)
    if n == 0:
        return 0
    if not triples:
        return n  # every element its own color is rainbow-free
    hi = n if upper is None else min(upper, n)
    for r in range(hi, 0, -1):
        if _exists_rainbow_free_with_r_colors_sat(n, triples, r):
            return r
    return 1


def rb_sat(orders: Sequence[int], include_zero: bool = True) -> int:
    """rb(G) via the SAT solver (scales past brute force)."""
    elems, triples = schur_triples(orders, include_zero)
    n = len(elems)
    if n == 0:
        return 1
    if not triples:
        return n + 1
    return max_rainbow_free_colors_sat(orders, include_zero) + 1


# --------------------------------------------------------------------------
# Convenience: dispatch (brute force for tiny, SAT otherwise) + group naming
# --------------------------------------------------------------------------

def rb(orders: Sequence[int], include_zero: bool = True,
       method: str = "auto") -> int:
    elems, _ = schur_triples(orders, include_zero)
    n = len(elems)
    if method == "brute":
        return rb_bruteforce(orders, include_zero)
    if method == "sat":
        return rb_sat(orders, include_zero)
    # auto: brute force is reliable and fast for n <= 11  (escape note: the
    # group-product symbol in docstrings is 'x', not a backslash sequence)
    if n <= 11:
        return rb_bruteforce(orders, include_zero)
    return rb_sat(orders, include_zero)


def group_name(orders: Sequence[int]) -> str:
    """Human-readable name like 'Z2 x Z2' (raw factor form, not invariant)."""
    return " x ".join(f"Z{n}" for n in orders)


# --------------------------------------------------------------------------
# Enumerating abelian groups up to isomorphism (primary decomposition)
# --------------------------------------------------------------------------

def abelian_groups(order: int) -> List[Tuple[int, ...]]:
    """All abelian groups of a given order, up to isomorphism, returned as
    sorted tuples of prime-power cyclic factors (the primary decomposition).
    Each tuple appears exactly once -> no isomorphic duplicates."""
    from sympy import factorint
    from sympy.utilities.iterables import partitions
    f = factorint(order)
    per_prime: List[List[List[int]]] = []
    for p, a in f.items():
        plist: List[List[int]] = []
        for part in partitions(a):
            powers: List[int] = []
            for size, mult in part.items():
                powers += [p ** size] * mult
            plist.append(sorted(powers))
        per_prime.append(plist)
    groups: List[Tuple[int, ...]] = []
    for combo in itertools.product(*per_prime):
        factors: List[int] = []
        for lst in combo:
            factors += lst
        groups.append(tuple(sorted(factors)))
    return groups


def group_rank(orders: Sequence[int]) -> int:
    """Rank = minimum number of generators = number of invariant factors.
    Computed from the primary decomposition: for each prime, the number of
    cyclic factors of that prime; the rank is the max over primes."""
    from sympy import factorint
    from collections import Counter
    by_prime: Counter = Counter()
    for q in orders:
        f = factorint(q)
        assert len(f) == 1, "factor must be a prime power"
        (p,) = f.keys()
        by_prime[p] += 1
    return max(by_prime.values()) if by_prime else 0


def group_exponent(orders: Sequence[int]) -> int:
    """Exponent = lcm of the orders of all elements = lcm of the factors."""
    import math
    e = 1
    for q in orders:
        e = e * q // math.gcd(e, q)
    return e


def is_cyclic(orders: Sequence[int]) -> bool:
    return group_rank(orders) <= 1


if __name__ == "__main__":
    import sys
    # quick smoke: cyclic small values, include_zero=True
    print("rb(Z_n, include_zero=True) for n=2..11 (brute force):")
    for n in range(2, 12):
        print(f"  Z{n}: {rb_bruteforce([n], include_zero=True)}")
