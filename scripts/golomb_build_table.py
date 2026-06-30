#!/usr/bin/env python3
"""Incrementally compute the G(m,n) table, writing each entry as it finishes.

Writes results to results_path as JSONL so partial progress survives a
timeout.  Cross-checks small cells against brute force.

Usage: python3 scripts/golomb_build_table.py MAXN [results.jsonl]
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.golomb_rectangles import (  # noqa: E402
    brute_force_max,
    golomb_rectangle,
    verify_witness,
)


def main() -> None:
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "scripts" / "golomb_results.jsonl"
    done = {}
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done[(r["m"], r["n"])] = r
    f = out.open("a")
    t_all = time.time()
    for m in range(1, maxn + 1):
        for n in range(m, maxn + 1):
            if (m, n) in done:
                print(f"  G({m},{n})={done[(m,n)]['G']} (cached)", flush=True)
                continue
            t0 = time.time()
            res = golomb_rectangle(m, n, symmetry=True)
            dt = time.time() - t0
            assert verify_witness(m, n, res["witness"]), f"bad witness {m},{n}"
            rec = {"m": m, "n": n, "G": res["G"], "proven": res["proven"],
                   "witness": res["witness"], "seconds": round(dt, 3)}
            # Brute-force cross-check where feasible.
            if m * n <= 28:
                bg, _ = brute_force_max(m, n)
                rec["brute"] = bg
                assert bg == res["G"], f"MISMATCH ({m},{n}): brute {bg} != sat {res['G']}"
                rec["brute_ok"] = True
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print(f"  G({m},{n})={res['G']:2d}  proven={res['proven']}  "
                  f"{dt:7.2f}s  brute={rec.get('brute','-')}", flush=True)
    f.close()
    print(f"TOTAL {time.time() - t_all:.1f}s", flush=True)


if __name__ == "__main__":
    main()
