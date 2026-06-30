"""Regenerate README.md tables from contributions.yaml (single source of truth).

Usage: python3 scripts/render_readme.py
Rewrites the region between the AUTOGEN markers in README.md.
"""
from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
BEGIN = "<!-- AUTOGEN:contributions BEGIN -->"
END = "<!-- AUTOGEN:contributions END -->"
OEIS = "https://oeis.org"


def link_id(sid: str) -> str:
    return f"[{sid}]({OEIS}/{sid})" if sid.startswith("A") and sid[1:].isdigit() else sid


def src_link(doc, s) -> str:
    base = doc["source_repos"].get(s["source_repo"], "")
    p = s.get("source_path", "")
    return f"[{s['source_repo']}/{p}]({base}/blob/main/{p})" if base and p else s.get("source_repo", "")


def table(doc, rows) -> str:
    out = ["| Sequence | Description | Provenance | Status |",
           "|----------|-------------|------------|--------|"]
    for s in rows:
        terms = (" " + ", ".join(s["my_terms"])) if s.get("my_terms") else ""
        out.append(f"| {link_id(s['id'])} | {s['title']}{terms} | {src_link(doc, s)} | {s['status']} |")
    return "\n".join(out)


def main() -> None:
    doc = yaml.safe_load((ROOT / "contributions.yaml").read_text())
    seqs = doc["sequences"]
    authored = [s for s in seqs if s["role"] == "authored" and s["status"] == "published"]
    extended = [s for s in seqs if s["role"] == "extended" and s["status"] == "published"]
    candidate = [s for s in seqs if s["status"] in ("candidate", "in-review")]
    other = [s for s in seqs if s["status"] in ("withdrawn", "not-submittable")]

    body = "\n\n".join([
        "### Authored (new sequences)\n\n" + table(doc, authored),
        "### Extended (new terms on classical sequences)\n\n" + table(doc, extended),
        "### Candidate / in review\n\n" + table(doc, candidate),
        "### Other dispositions\n\n" + table(doc, other),
    ])

    readme = (ROOT / "README.md").read_text()
    pre = readme.split(BEGIN)[0]
    post = readme.split(END)[1] if END in readme else "\n"
    (ROOT / "README.md").write_text(f"{pre}{BEGIN}\n\n{body}\n\n{END}{post}")
    print(f"README.md regenerated: {len(authored)} authored, {len(extended)} extended, "
          f"{len(candidate)} candidate/in-review.")


if __name__ == "__main__":
    main()
