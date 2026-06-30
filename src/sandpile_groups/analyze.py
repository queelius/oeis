"""Analyze the sandpile-group catalog for patterns and OEIS matches.

Reports:
  - Sequence of distinct sandpile groups per n (search OEIS)
  - Sequence of cyclic sandpile groups per n (search OEIS)
  - Maximum group order per n (= n^{n-2}, A000272)
  - Trees count (A000055)
  - For each small group order k, count graphs with |K(G)|=k by n
  - For each small group structure, count graphs with K(G) = G_target by n
  - Rank distribution
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

OUT_DIR = Path(__file__).parent
CATALOG = OUT_DIR / "catalog.jsonl"


def load() -> list[dict]:
    out = []
    with CATALOG.open() as f:
        for line in f:
            out.append(json.loads(line))
    return out


def main():
    rows = load()
    by_n: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        by_n[r["n"]].append(r)
    ns = sorted(by_n)

    print(f"Loaded {len(rows)} records across n={ns[0]}..{ns[-1]}")
    print()

    # 1. Counts per n
    print("=== Per-n summary ===")
    print("n  | graphs | distinct K | trees | cyclic | max|K| | min|K|>1")
    print("---+--------+------------+-------+--------+--------+---------")
    for n in ns:
        g = by_n[n]
        sigs = {r["signature"] for r in g}
        trees = sum(1 for r in g if r["is_tree"])
        cyclic = sum(1 for r in g if r["is_cyclic"])
        max_o = max(r["order"] for r in g)
        nontriv = [r["order"] for r in g if r["order"] > 1]
        min_o = min(nontriv) if nontriv else None
        print(f"{n}  | {len(g):6d} | {len(sigs):10d} | {trees:5d} | "
              f"{cyclic:6d} | {max_o:6d} | {min_o}")
    print()

    # 2. OEIS-hunt sequences
    print("=== Candidate OEIS sequences ===")
    seq_distinct = [len({r["signature"] for r in by_n[n]}) for n in ns]
    seq_cyclic = [sum(1 for r in by_n[n] if r["is_cyclic"]) for n in ns]
    seq_trees = [sum(1 for r in by_n[n] if r["is_tree"]) for n in ns]
    seq_max_order = [max(r["order"] for r in by_n[n]) for n in ns]
    seq_total = [len(by_n[n]) for n in ns]
    seq_noncyclic = [
        sum(1 for r in by_n[n] if not r["is_cyclic"])
        for n in ns
    ]
    seq_unicyclic = [
        sum(1 for r in by_n[n] if r["m"] == r["n"])  # connected, m=n => unicyclic
        for n in ns
    ]
    print(f"  ns                                = {ns}")
    print(f"  distinct sandpile groups          = {seq_distinct}")
    print(f"    -> OEIS check (search exact)")
    print(f"  cyclic sandpile groups            = {seq_cyclic}")
    print(f"    -> OEIS check (search exact)")
    print(f"  non-cyclic sandpile groups        = {seq_noncyclic}")
    print(f"  trees (= |K|=1)                   = {seq_trees}")
    print(f"    -> A000055 expected: 1,1,1,2,3,6,11,23")
    print(f"  max |K| (= n^{{n-2}})            = {seq_max_order}")
    print(f"    -> A000272 expected: 1,1,3,16,125,1296,16807,262144")
    print(f"  total connected graphs            = {seq_total}")
    print(f"    -> A001349 expected: 1,1,2,6,21,112,853,11117")
    print(f"  unicyclic connected graphs        = {seq_unicyclic}")
    print(f"    -> A001429 expected: 0,0,1,2,5,13,33,89")
    print()

    # 3. Counts of K(G) = Z/k for small k.
    print("=== Counts of cyclic K(G) = Z/k by n (small k) ===")
    for k in range(1, 13):
        row = []
        for n in ns:
            cnt = sum(1 for r in by_n[n] if r["is_cyclic"] and r["order"] == k)
            row.append(cnt)
        print(f"  K = Z/{k:<3d}: {row}")
    print()

    # 4. Most common sandpile groups overall
    print("=== Top 25 most common sandpile groups across catalog ===")
    sig_counts: Counter[str] = Counter(r["signature"] for r in rows)
    for sig, cnt in sig_counts.most_common(25):
        # Per-n breakdown
        by_n_cnts = Counter(r["n"] for r in rows if r["signature"] == sig)
        n_breakdown = ", ".join(f"n={n}:{by_n_cnts[n]}" for n in ns if by_n_cnts[n])
        print(f"  {sig:30s} total={cnt:5d}  ({n_breakdown})")
    print()

    # 5. Rank distribution per n
    print("=== Rank distribution per n (rank = # nontrivial invariant factors) ===")
    print("n  | rank=0 (tree) | rank=1 (cyclic) | rank=2 | rank=3 | rank=4 | rank=5 | rank=6")
    for n in ns:
        ranks = Counter(r["rank"] for r in by_n[n])
        cells = [str(ranks.get(k, 0)) for k in range(7)]
        print(f"{n}  | " + " | ".join(c.rjust(13 if i==0 else (15 if i==1 else 6)) for i, c in enumerate(cells)))
    print()

    # 6. Maximum rank per n -- when is K(G) achieved at full rank n-1?
    # Full-rank invariant factor sequence means d_1 > 1 -- the group is "highly
    # non-cyclic".  K_n achieves rank n-2 (since (Z/n)^{n-2}).
    print("=== Maximum rank per n ===")
    for n in ns:
        max_rank = max(r["rank"] for r in by_n[n])
        # Which graphs?
        achievers = [r for r in by_n[n] if r["rank"] == max_rank]
        # Identify K_n
        kn_edges = n * (n - 1) // 2
        is_kn = any(r["m"] == kn_edges for r in achievers)
        print(f"  n={n}: max rank={max_rank}, achievers={len(achievers)} "
              f"(K_n among them: {is_kn})")
    print()

    # 7. Number of distinct group orders per n
    print("=== Distinct |K(G)| values per n ===")
    for n in ns:
        orders = {r["order"] for r in by_n[n]}
        print(f"  n={n}: {len(orders)} distinct orders, range [{min(orders)}, {max(orders)}]")
    print()

    # 8. Count of graphs with K(G) trivial (= trees) vs cyclic non-trivial
    #    vs non-cyclic, per n
    print("=== Tree / cyclic / non-cyclic per n ===")
    print("n  | trees | cyclic-nontrivial | non-cyclic")
    for n in ns:
        g = by_n[n]
        trees = sum(1 for r in g if r["is_tree"])
        cyc_nt = sum(1 for r in g if r["is_cyclic"] and not r["is_tree"])
        non_cyc = sum(1 for r in g if not r["is_cyclic"])
        print(f"{n}  | {trees:5d} | {cyc_nt:18d} | {non_cyc:10d}")
    print()


if __name__ == "__main__":
    main()
