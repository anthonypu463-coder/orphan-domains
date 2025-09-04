#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

IN_CONF = Path("data/outputs/hmmscan/pfam_hits_confirmed.tsv")
OUT_KEEP = Path("data/outputs/hmmscan/pfam_domains_nonoverlap.tsv")
OUT_DROP = Path("data/outputs/hmmscan/pfam_domains_dropped.tsv")
PFAM_CLANS_TSV = Path("data/pfam/Pfam-A.clans.tsv")  # optional; if absent, clan precedence is ignored

def norm_acc(acc: str) -> str:
    # Normalize PFxxxxx[.yy] → PFxxxxx
    if not acc:
        return ""
    i = acc.find(".")
    return acc[:i] if i > 0 else acc

def load_clans(path: Path) -> Dict[str, str]:
    """
    Load mapping PFxxxxx(.yy) -> CLxxxxx (if file present). Tolerant to unknown/extra columns.
    Returns empty dict if file missing.
    """
    if not path.exists():
        return {}
    m: Dict[str, str] = {}
    with path.open() as fh:
        for raw in fh:
            ln = raw.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = [p for p in ln.replace(",", "\t").split() if p]  # handle tab- or comma-delimited
            pfam = next((p for p in parts if p.startswith("PF")), "")
            clan = next((p for p in parts if p.startswith("CL")), "")
            if pfam and clan:
                m[norm_acc(pfam)] = clan
    return m

def overlap_len(a1: int, a2: int, b1: int, b2: int) -> int:
    lo = max(a1, b1); hi = min(a2, b2)
    return max(0, hi - lo + 1)

def parse_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except Exception:
        return None

def read_hits(infile: Path) -> List[dict]:
    if not infile.exists():
        raise SystemExit(f"Missing {infile}")
    recs: List[dict] = []
    with infile.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            # Required fields; robust to blanks
            recs.append({
                "uniprot": row["uniprot"],
                "pfam_id": row["pfam_id"],
                "pfam_acc": row["pfam_acc"],
                "qlen": int(row["qlen"]) if row["qlen"].isdigit() else 0,
                "env_from": int(row["env_from"]) if row["env_from"].isdigit() else 0,
                "env_to": int(row["env_to"]) if row["env_to"].isdigit() else 0,
                "dom_score": parse_float(row["dom_score"]),
                "full_score": parse_float(row["full_score"]),
                "description": row.get("description", ""),
            })
    return recs

def resolve_for_protein(recs: List[dict], clans: Dict[str, str], trim_small: bool = False) -> Tuple[List[dict], List[dict]]:
    """
    Greedy, score-descending selection:
    - Sort by dom_score desc, then full_score desc, then longer domain first.
    - Accept first interval; for each next one, drop if overlaps any accepted interval.
    - If trim_small=True and overlap <= 10 residues or <=10% of smaller domain, trim loser.
    """
    # Sort by significance
    recs = sorted(
        recs,
        key=lambda r: (
            -(r["dom_score"] if r["dom_score"] is not None else float("-inf")),
            -(r["full_score"] if r["full_score"] is not None else float("-inf")),
            -(r["env_to"] - r["env_from"] + 1),
        ),
    )
    kept: List[dict] = []
    dropped: List[dict] = []

    for r in recs:
        cur_l, cur_r = r["env_from"], r["env_to"]
        if cur_l <= 0 or cur_r <= 0 or cur_l > cur_r:
            r2 = {**r, "reason": "invalid_coords"}
            dropped.append(r2)
            continue

        overlaps = [k for k in kept if overlap_len(cur_l, cur_r, k["env_from"], k["env_to"]) > 0]
        if not overlaps:
            kept.append(r)
            continue

        # With greedy score-sorted order, any overlap implies this record scores <= existing ones.
        # Clan precedence (if same clan): keep higher-scoring (the existing kept one by construction).
        # Different clans: default policy = drop the later, optional trim for tiny overlaps.
        if not trim_small:
            dropped.append({**r, "reason": "overlap_drop"})
            continue

        # Attempt trimming for tiny overlaps
        # Compute total overlap with all kept regions; if resulting remaining span is meaningful, keep trimmed.
        covered = []
        for k in overlaps:
            covered.append((max(cur_l, k["env_from"]), min(cur_r, k["env_to"])))
        # Merge covered intervals
        covered.sort()
        merged = []
        for s, e in covered:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        # Subtract merged from current interval -> at most two segments remain
        segments: List[Tuple[int,int]] = []
        start = cur_l
        for s,e in merged:
            if start < s:
                segments.append((start, s-1))
            start = e + 1
        if start <= cur_r:
            segments.append((start, cur_r))

        # Keep the longest remaining segment if it's non-trivial
        if segments:
            seg = max(segments, key=lambda x: x[1]-x[0])
            seg_len = seg[1] - seg[0] + 1
            base_len = cur_r - cur_l + 1
            tiny = seg_len <= 10 or seg_len / base_len <= 0.10
            if not tiny:
                kept.append({**r, "env_from": seg[0], "env_to": seg[1], "note": "trimmed"})
                continue

        dropped.append({**r, "reason": "overlap_drop"})
    return kept, dropped

def main() -> None:
    all_hits = read_hits(IN_CONF)
    clans = load_clans(PFAM_CLANS_TSV)
    by_prot: Dict[str, List[dict]] = {}
    for r in all_hits:
        by_prot.setdefault(r["uniprot"], []).append(r)

    kept_all: List[dict] = []
    drop_all: List[dict] = []
    for acc, lst in by_prot.items():
        kept, drop = resolve_for_protein(lst, clans, trim_small=False)
        kept_all.extend(kept)
        for d in drop:
            d["uniprot"] = acc  # ensure present
        drop_all.extend(drop)

    # Write outputs
    OUT_KEEP.parent.mkdir(parents=True, exist_ok=True)
    with OUT_KEEP.open("w") as fh:
        fh.write("\t".join(["uniprot","pfam_id","pfam_acc","env_from","env_to","dom_score","full_score","description"]) + "\n")
        for r in kept_all:
            fh.write("\t".join([
                r["uniprot"], r["pfam_id"], r["pfam_acc"],
                str(r["env_from"]), str(r["env_to"]),
                "" if r["dom_score"] is None else f"{r['dom_score']:.2f}",
                "" if r["full_score"] is None else f"{r['full_score']:.2f}",
                r.get("description",""),
            ]) + "\n")

    with OUT_DROP.open("w") as fh:
        fh.write("\t".join(["uniprot","pfam_id","pfam_acc","env_from","env_to","dom_score","full_score","reason"]) + "\n")
        for r in drop_all:
            fh.write("\t".join([
                r["uniprot"], r["pfam_id"], r["pfam_acc"],
                str(r["env_from"]), str(r["env_to"]),
                "" if r["dom_score"] is None else f"{r['dom_score']:.2f}",
                "" if r["full_score"] is None else f"{r['full_score']:.2f}",
                r.get("reason",""),
            ]) + "\n")

    print(f"Wrote {OUT_KEEP} ({sum(1 for _ in OUT_KEEP.open())-1} rows)")
    print(f"Wrote {OUT_DROP} ({sum(1 for _ in OUT_DROP.open())-1} rows)")

if __name__ == "__main__":
    main()
