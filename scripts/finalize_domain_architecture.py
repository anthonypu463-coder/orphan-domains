#!/usr/bin/env python3
from __future__ import annotations
import csv, argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

CONF = Path("data/outputs/hmmscan/pfam_hits_confirmed.tsv")
MANI = Path("data/afdb/index/reference_manifest.tsv")
CLANS = Path("data/pfam/Pfam-A.clans.tsv")   # optional

OUT_DIR = Path("data/outputs/architecture")
OUT_DOM = OUT_DIR / "pfam_domains_final.tsv"
OUT_ARCH= OUT_DIR / "pfam_architecture_final.tsv"

def norm_acc(acc: str) -> str:
    if not acc: return ""
    i = acc.find(".")
    return acc[:i] if i > 0 else acc

def load_clans(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    m: Dict[str, str] = {}
    with path.open() as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            # Pfam-A.clans.tsv is tab-separated: <pfam_acc> <clan_acc> <pfam_id> <clan_id> ...
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].startswith("PF") and parts[1].startswith("CL"):
                m[norm_acc(parts[0])] = parts[1]
    return m

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0

def to_float(s: str) -> Optional[float]:
    try: return float(s)
    except Exception: return None

def read_confirmed(path: Path) -> List[dict]:
    if not path.exists():
        raise SystemExit(f"Missing {path}")
    recs: List[dict] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            recs.append({
                "uniprot": r["uniprot"],
                "pfam_id": r["pfam_id"],
                "pfam_acc": r["pfam_acc"],
                "env_from": to_int(r["env_from"]),
                "env_to": to_int(r["env_to"]),
                "dom_score": to_float(r["dom_score"]),
                "full_score": to_float(r["full_score"]),
                "description": r.get("description",""),
            })
    return recs

def read_lengths(path: Path) -> Dict[str, int]:
    if not path.exists():
        raise SystemExit(f"Missing {path}")
    m: Dict[str,int] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r.get("uniprot","").strip()
            if acc:
                m[acc] = to_int(r.get("residues","0"))
    return m

def overlap(a1:int,a2:int,b1:int,b2:int) -> int:
    lo, hi = max(a1,b1), min(a2,b2)
    return max(0, hi-lo+1)

def covered_length(intervals: List[Tuple[int,int]]) -> int:
    if not intervals: return 0
    iv = sorted(intervals)
    total = 0
    cur_l, cur_r = iv[0]
    for s,e in iv[1:]:
        if s <= cur_r:
            cur_r = max(cur_r, e)
        else:
            total += (cur_r - cur_l + 1)
            cur_l, cur_r = s, e
    total += (cur_r - cur_l + 1)
    return total

def resolve_for_protein(
    recs: List[dict],
    clans: Dict[str,str],
    trim_policy: str = "same-clan",
    min_fragment: int = 30
) -> List[dict]:
    """
    Greedy keep of highest-scoring domains; optional trimming of small overlaps.
    trim_policy: "none" | "same-clan" | "any"
    """
    recs = [r for r in recs if r["env_from"]>0 and r["env_to"]>=r["env_from"]]
    # Score sort: dom_score desc, full_score desc, length desc
    recs.sort(key=lambda r: (
        -(r["dom_score"] if r["dom_score"] is not None else float("-inf")),
        -(r["full_score"] if r["full_score"] is not None else float("-inf")),
        -(r["env_to"] - r["env_from"] + 1),
    ))
    kept: List[dict] = []
    for r in recs:
        l, rgt = r["env_from"], r["env_to"]
        overlaps = [k for k in kept if overlap(l, rgt, k["env_from"], k["env_to"]) > 0]
        if not overlaps:
            kept.append(r); continue

        # Determine if trimming is permitted
        allow_trim = False
        if trim_policy == "any":
            allow_trim = True
        elif trim_policy == "same-clan":
            clan_r = clans.get(norm_acc(r["pfam_acc"]), "")
            # permitted if every overlap is same clan
            allow_trim = bool(clan_r) and all(clans.get(norm_acc(k["pfam_acc"]), "") == clan_r for k in overlaps)

        if not allow_trim:
            # lower-ranked by construction; drop it
            continue

        # Compute parts of r not covered by kept
        cov = []
        for k in overlaps:
            cov.append((max(l, k["env_from"]), min(rgt, k["env_to"])))
        cov.sort()
        merged = []
        for s,e in cov:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s,e])
            else:
                merged[-1][1] = max(merged[-1][1], e)

        segments: List[Tuple[int,int]] = []
        cursor = l
        for s,e in merged:
            if cursor < s:
                segments.append((cursor, s-1))
            cursor = e+1
        if cursor <= rgt:
            segments.append((cursor, rgt))

        if not segments:
            continue  # fully covered by better domains

        # Keep the longest remaining fragment if it is meaningful
        seg = max(segments, key=lambda x: x[1]-x[0])
        if (seg[1]-seg[0]+1) >= min_fragment:
            kept.append({**r, "env_from": seg[0], "env_to": seg[1], "note": "trimmed"})
        # else: drop tiny remainder
    # Sort kept by coordinate for readability
    kept.sort(key=lambda x: (x["env_from"], x["env_to"]))
    return kept

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trim-policy", choices=["none","same-clan","any"], default="same-clan",
                    help="Allow trimming only for same-clan overlaps (default), for any overlaps, or disable.")
    ap.add_argument("--min-fragment", type=int, default=30, help="Minimum length to keep a trimmed fragment.")
    args = ap.parse_args()

    clans = load_clans(CLANS)
    confirmed = read_confirmed(CONF)
    lengths = read_lengths(MANI)

    by = {}
    for r in confirmed:
        by.setdefault(r["uniprot"], []).append(r)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Write per-domain final map
    with OUT_DOM.open("w") as fh:
        fh.write("\t".join(["uniprot","pfam_id","pfam_acc","clan","start","end","dom_score","full_score","note"]) + "\n")
        for acc, lst in by.items():
            kept = resolve_for_protein(lst, clans, trim_policy=args.trim_policy, min_fragment=args.min_fragment)
            for k in kept:
                clan = clans.get(norm_acc(k["pfam_acc"]), "")
                fh.write("\t".join([
                    acc, k["pfam_id"], k["pfam_acc"], clan,
                    str(k["env_from"]), str(k["env_to"]),
                    "" if k["dom_score"] is None else f"{k['dom_score']:.2f}",
                    "" if k["full_score"] is None else f"{k['full_score']:.2f}",
                    k.get("note",""),
                ]) + "\n")

    # Summarize per-protein architecture and coverage
    # Re-read OUT_DOM to avoid keeping all in memory
    spans: Dict[str, List[Tuple[int,int]]] = {}
    dom_counts: Dict[str,int] = {}
    with OUT_DOM.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            s, e = to_int(r["start"]), to_int(r["end"])
            if s>0 and e>=s:
                spans.setdefault(acc, []).append((s,e))
            dom_counts[acc] = dom_counts.get(acc, 0) + 1

    with OUT_ARCH.open("w") as fh:
        fh.write("\t".join(["uniprot","residues","n_domains","covered_residues","coverage_frac"]) + "\n")
        for acc, L in lengths.items():
            iv = spans.get(acc, [])
            cov = covered_length(iv)
            frac = (cov / L) if L>0 else 0.0
            fh.write(f"{acc}\t{L}\t{dom_counts.get(acc,0)}\t{cov}\t{frac:.4f}\n")

    print(f"Wrote {OUT_DOM}")
    print(f"Wrote {OUT_ARCH}")

if __name__ == "__main__":
    main()
