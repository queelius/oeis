#!/usr/bin/env python3
"""Reproduce / verify OEIS A395644 from the committed KnotInfo-derived data.

A395644: a(n) = number of fibered prime knots with n crossings, n >= 3.
Offset 3,1. Published DATA (11 terms, n = 3..13):

    1, 1, 1, 2, 3, 12, 23, 74, 256, 873, 4151

Data provenance and scope
-------------------------
The published sequence was authored from the canonical KnotInfo database
(Livingston-Moore), counting prime knots (up to mirror image, Perko pair
excluded) whose `fibered` flag is set, by crossing number, for n = 3..13.

The committed snapshot data/catalog.json (801 prime knots) was generated in
the source repo `open-problems` by src/knot_invariants.py via the
`spherogram` + `knot_floer_homology` packages, and only covers n <= 11.
It therefore can speak to a(3)..a(11) but NOT a(12), a(13): the n=12 (2176
knots) and n=13 (9988 knots) fibered counts require the KnotInfo database
directly and are outside the committed catalog. This is the KnotInfo data
dependency for the tail of the sequence.

The n=11 anomaly (handled below)
--------------------------------
For n = 3..10 the catalog reproduces the published values exactly (this
region passed 105/105 KnotInfo per-knot cross-checks in the source repo).
At n = 11 the raw catalog gives 257 fibered knots, one MORE than the
published 256. The source repo's c11-anomaly-investigation.md diagnosed
this precisely: spherogram's c=11 table includes BOTH chiralities of the
knot KnotInfo records once as 11n_128 (our names 11_376 and 11_449, mirror
images, both fibered), while omitting one chirality of another knot. The
total stays 552 = A002863(11) by coincidental cancellation, but the
fibered column is inflated by exactly one. Counting up to mirror image (the
KnotInfo / OEIS convention) removes this one duplicated fibered chiral
knot, recovering a(11) = 256.

This script is SELF-CONTAINED: it does NOT import spherogram or
database_knotinfo. It recomputes a(3)..a(11) from the committed catalog,
applies the documented mirror-deduplication at n=11, asserts an exact match
to the published prefix, and reports the KnotInfo dependency for a(12),
a(13).

Runtime: well under 1 second.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

# Full published DATA (n = 3..13).
PUBLISHED = [1, 1, 1, 2, 3, 12, 23, 74, 256, 873, 4151]
PUBLISHED_N = list(range(3, 3 + len(PUBLISHED)))  # 3..13

# spherogram's c=11 table double-counts this mirror pair (KnotInfo: 11n_128,
# listed once). Both are fibered; counting up to mirror image removes one.
SPHEROGRAM_C11_CHIRAL_DUP = {"11_376", "11_449"}

CATALOG = Path(__file__).parent / "data" / "catalog.json"


def recompute_prefix() -> dict[int, int]:
    """Fibered prime-knot count per crossing number n=3..11 from the catalog,
    with the documented n=11 mirror-deduplication applied (KnotInfo convention).
    """
    catalog = json.loads(CATALOG.read_text())["catalog"]

    fibered = Counter(d["crossings"] for d in catalog if d["fibered"])

    # Apply mirror-dedup at n=11: spherogram lists both chiralities of the
    # single KnotInfo knot 11n_128 (both fibered) -> subtract one.
    present = {d["name"] for d in catalog}
    if SPHEROGRAM_C11_CHIRAL_DUP <= present:
        fibered[11] -= 1

    return {n: fibered.get(n, 0) for n in range(3, 12)}


def main() -> None:
    recomputed = recompute_prefix()

    # Compare against the portion of PUBLISHED the catalog (n<=11) can reach.
    reachable_n = [n for n in PUBLISHED_N if n <= 11]
    for n in reachable_n:
        got = recomputed[n]
        want = PUBLISHED[PUBLISHED_N.index(n)]
        assert got == want, (
            f"n={n}: recomputed fibered count {got} != published {want}"
        )

    # Sanity: total prime knots per crossing must match A002863 (n<=11).
    catalog = json.loads(CATALOG.read_text())["catalog"]
    a002863 = {3: 1, 4: 1, 5: 2, 6: 3, 7: 7, 8: 21, 9: 49, 10: 165, 11: 552}
    totals = Counter(d["crossings"] for d in catalog)
    for n, expected in a002863.items():
        assert totals[n] == expected, f"n={n} total {totals[n]} != A002863 {expected}"

    tail_n = [n for n in PUBLISHED_N if n > 11]  # 12, 13
    print(
        f"OK A395644: recomputed a(3..11) = "
        f"{[recomputed[n] for n in range(3, 12)]} from catalog.json "
        f"(mirror-dedup applied at n=11: spherogram double-counts chiral "
        f"11n_128); exact match to published prefix. "
        f"Tail a(12)={PUBLISHED[PUBLISHED_N.index(12)]}, "
        f"a(13)={PUBLISHED[PUBLISHED_N.index(13)]} require the KnotInfo "
        f"database (n={tail_n} beyond the c<=11 committed catalog)."
    )


if __name__ == "__main__":
    main()
