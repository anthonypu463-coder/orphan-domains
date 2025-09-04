#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import re, sys

"""
Build a GA threshold table from Pfam-A.hmm (HMMER3 text format).
Outputs: data/outputs/hmmscan/pfam_ga.tsv with columns:
pfam_acc, pfam_id, ga_seq, ga_dom
Notes:
- Pfam uses NAME (not ID) in HMMER3 headers; we also accept ID if present.
- GA lines look like: "GA 25.0 25.0;"
"""

def main() -> None:
    hmm = Path("data/pfam/Pfam-A.hmm")
    out = Path("data/outputs/hmmscan/pfam_ga.tsv")
    if not hmm.exists():
        sys.exit(f"Missing {hmm}")

    name = acc = None
    ga_seq = ga_dom = None
    rows = []

    re_acc  = re.compile(r"^ACC\s+(\S+)")
    re_name = re.compile(r"^NAME\s+(\S+)")
    re_id   = re.compile(r"^ID\s+(\S+)")
    re_ga   = re.compile(r"^GA\s+([0-9.]+)\s+([0-9.]+);")

    with hmm.open() as fh:
        for ln in fh:
            if ln.startswith("ACC"):
                m = re_acc.match(ln); acc = m.group(1) if m else acc
            elif ln.startswith("NAME"):
                m = re_name.match(ln); name = m.group(1) if m else name
            elif ln.startswith("ID"):
                # Some files may still carry ID; prefer NAME, but fall back if NAME absent.
                if not name:
                    m = re_id.match(ln); name = m.group(1) if m else name
            elif ln.startswith("GA"):
                m = re_ga.match(ln)
                if m:
                    ga_seq = float(m.group(1)); ga_dom = float(m.group(2))
            elif ln.strip() == "//":
                if acc and name and ga_seq is not None and ga_dom is not None:
                    rows.append((acc, name, ga_seq, ga_dom))
                name = acc = None
                ga_seq = ga_dom = None

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("pfam_acc\tpfam_id\tga_seq\tga_dom\n")
        for a, n, gs, gd in rows:
            fh.write(f"{a}\t{n}\t{gs}\t{gd}\n")
    print(f"Wrote {out} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
