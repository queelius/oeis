"""Build the sandpile-group catalog for connected graphs on n <= 8.

Outputs:
  catalog.jsonl      one JSON record per connected graph
  summary.json       per-n counts, distinct groups, OEIS-friendly stats

Run from the repository root:

    python3 discoveries/algebraic-graph-theory/sandpile-groups/build_catalog.py
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

# Make the project src/ importable when run from anywhere.
_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from src.sandpile import (  # noqa: E402
    catalog_record,
    enumerate_connected_graphs,
    is_cyclic,
    num_spanning_trees,
)


OUT_DIR = Path(__file__).parent
CATALOG = OUT_DIR / "catalog.jsonl"
SUMMARY = OUT_DIR / "summary.json"


def main(n_max: int = 8) -> None:
    summary = {
        "n_max": n_max,
        "by_n": {},
    }
    with CATALOG.open("w") as fout:
        total = 0
        for n in range(2, n_max + 1):
            t0 = time.time()
            count = 0
            distinct_groups: Counter[str] = Counter()
            order_counts: Counter[int] = Counter()
            cyclic = 0
            tree = 0
            ranks: Counter[int] = Counter()
            max_order = 0
            min_nontrivial_order = None
            for (_, edges) in enumerate_connected_graphs(n):
                rec = catalog_record(n, edges)
                # Cross-check |K(G)| against matrix-tree (only for small n
                # to keep things fast).
                if n <= 6:
                    nst = num_spanning_trees(
                        [[1 if (i, j) in {tuple(sorted(e)) for e in edges}
                          or (j, i) in {tuple(sorted(e)) for e in edges} else 0
                          for j in range(n)] for i in range(n)]
                    )
                    assert nst == rec["order"], (n, edges, nst, rec)
                fout.write(json.dumps(rec) + "\n")
                count += 1
                total += 1
                distinct_groups[rec["signature"]] += 1
                order_counts[rec["order"]] += 1
                ranks[rec["rank"]] += 1
                if rec["is_cyclic"]:
                    cyclic += 1
                if rec["is_tree"]:
                    tree += 1
                if rec["order"] > max_order:
                    max_order = rec["order"]
                if rec["order"] > 1:
                    if min_nontrivial_order is None or rec["order"] < min_nontrivial_order:
                        min_nontrivial_order = rec["order"]
            dt = time.time() - t0
            summary["by_n"][str(n)] = {
                "graphs": count,
                "distinct_groups": len(distinct_groups),
                "trees": tree,
                "cyclic": cyclic,
                "max_order": max_order,
                "min_nontrivial_order": min_nontrivial_order,
                "rank_histogram": dict(sorted(ranks.items())),
                "top_orders": dict(order_counts.most_common(10)),
                "elapsed_seconds": round(dt, 3),
            }
            print(
                f"n={n}: {count:>6d} graphs, "
                f"{len(distinct_groups):>4d} distinct groups, "
                f"{tree:>4d} trees, {cyclic:>5d} cyclic, "
                f"max|K|={max_order} ({dt:.1f}s)"
            )
        summary["total_graphs"] = total
    with SUMMARY.open("w") as fout:
        json.dump(summary, fout, indent=2)
    print(f"\nWrote {CATALOG} ({total} records)")
    print(f"Wrote {SUMMARY}")


if __name__ == "__main__":
    n_max = 8
    if len(sys.argv) > 1:
        n_max = int(sys.argv[1])
    main(n_max)
