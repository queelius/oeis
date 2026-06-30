"""Lint contributions.yaml: required fields, A-number format, data_dir presence.

Usage: python3 scripts/validate.py
Exit 0 if clean, 1 otherwise.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
REQUIRED = {"id", "title", "role", "status", "domain", "source_repo"}
ROLES = {"authored", "extended", "comment"}
STATUSES = {"published", "in-review", "candidate", "withdrawn", "not-submittable"}
ANUM = re.compile(r"^A\d{6}$")


def main() -> int:
    doc = yaml.safe_load((ROOT / "contributions.yaml").read_text())
    seqs = doc.get("sequences", [])
    errors: list[str] = []
    ids = set()

    for s in seqs:
        sid = s.get("id", "<no id>")
        if sid in ids:
            errors.append(f"{sid}: duplicate id")
        ids.add(sid)
        missing = REQUIRED - set(s)
        if missing:
            errors.append(f"{sid}: missing fields {sorted(missing)}")
        if s.get("role") not in ROLES:
            errors.append(f"{sid}: bad role {s.get('role')!r}")
        if s.get("status") not in STATUSES:
            errors.append(f"{sid}: bad status {s.get('status')!r}")
        # published/in-review extensions and authored A-numbers must be A######
        if s.get("status") in ("published", "in-review") and sid.startswith("A"):
            if not ANUM.match(sid):
                errors.append(f"{sid}: malformed A-number")
        # data_dir, if given, must exist
        dd = s.get("data_dir")
        if dd and not (ROOT / dd).is_dir():
            errors.append(f"{sid}: data_dir {dd} missing")
        # bfile, if given, must exist
        bf = s.get("bfile")
        if bf and not (ROOT / bf).is_file():
            errors.append(f"{sid}: bfile {bf} missing")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print("  -", e)
        return 1
    print(f"OK: {len(seqs)} contribution records valid "
          f"({sum(1 for s in seqs if s['status']=='published')} published).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
