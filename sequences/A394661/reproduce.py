#!/usr/bin/env python3
"""Reproduce OEIS A394661 from the committed KnotInfo-derived catalog.

A394661: Triangle read by rows T(n, k) = number of prime knots with n
crossings and three-genus (Seifert genus) k, for n >= 3 and
1 <= k <= floor((n-1)/2). Offset 3,6. Perko pair (10_161 = 10_162)
excluded; knots counted up to mirror image, following KnotInfo.

Data provenance
---------------
The per-knot invariant catalog (data/catalog.json, 801 prime knots through
11 crossings) was generated in the source repo `open-problems` by
src/knot_invariants.py, which enumerates prime knots via the `spherogram`
package (Rolfsen table c <= 10, Hoste-Thistlethwaite c = 11) and computes
the Seifert genus per knot via knot_floer_homology (HFK). The genus column
was cross-checked against KnotInfo (Livingston-Moore) for all 35 prime
knots with c <= 8 (exact) and the per-crossing genus distribution matches
KnotInfo / the published sequence at every n = 3..11.

This script is SELF-CONTAINED: it does NOT import spherogram or
database_knotinfo. It recomputes the full triangle by re-tallying the
committed catalog.json (crossing, genus) and asserts an exact match
against the 25 published DATA terms.

Runtime: well under 1 second.
"""
from __future__ import annotations

import json
from collections import defaultdict, Counter
from pathlib import Path

# Published DATA (25 terms, rows n = 3..11), from the live OEIS entry.
PUBLISHED = [
    1,
    1,
    1, 1,
    1, 2,
    2, 4, 1,
    2, 10, 9,
    4, 22, 22, 1,
    2, 44, 93, 26,
    6, 96, 289, 160, 1,
]

CATALOG = Path(__file__).parent / "data" / "catalog.json"


def recompute_triangle() -> list[int]:
    """Re-tally the genus-by-crossing triangle from the committed catalog."""
    catalog = json.loads(CATALOG.read_text())["catalog"]

    genus_by_crossing: dict[int, Counter] = defaultdict(Counter)
    for knot in catalog:
        genus_by_crossing[knot["crossings"]][knot["genus"]] += 1

    flat: list[int] = []
    for n in sorted(genus_by_crossing):
        kmax = (n - 1) // 2  # max realized genus = floor((n-1)/2)
        flat.extend(genus_by_crossing[n].get(k, 0) for k in range(1, kmax + 1))
    return flat, genus_by_crossing


def main() -> None:
    flat, gbc = recompute_triangle()

    # Row sums must equal A002863 (prime knots by crossing number).
    a002863 = {3: 1, 4: 1, 5: 2, 6: 3, 7: 7, 8: 21, 9: 49, 10: 165, 11: 552}
    for n, expected in a002863.items():
        got = sum(gbc[n].values())
        assert got == expected, f"row n={n} sum {got} != A002863 {expected}"

    assert flat == PUBLISHED, (
        "recomputed triangle does not match published DATA:\n"
        f"  recomputed: {flat}\n"
        f"  published : {PUBLISHED}"
    )

    print(
        f"OK A394661: recomputed {len(flat)} triangle terms (rows n=3..11) "
        f"from {len(json.loads(CATALOG.read_text())['catalog'])} prime knots "
        f"in catalog.json; exact match to published DATA; row sums = A002863."
    )


if __name__ == "__main__":
    main()
