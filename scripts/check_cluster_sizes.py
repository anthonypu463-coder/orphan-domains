#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List

# Inputs from 5B/5C
REPS_TSV  = Path("data/outputs/cluster26k/cluster26k_representatives.tsv")  # cluster_id, representative, n_members, ...
MEMB_TSV  = Path("data/outputs/cluster26k/cluster26k_members.tsv")          # cluster_id, member
FASTA26K  = Path("data/outputs/orphans/orphan_strict.fasta")

# Outputs
OUT_DIR   = Path("data/outputs/cluster26k")
HIST_TSV  = OUT_DIR / "cluster_size_hist.tsv"
LARGE_TSV = OUT_DIR / "clusters_large.tsv"
SUMMARY   = OUT_DIR / "cluster_size_summary.json"

LARGE_THRESHOLD = 100  # flag clusters with size > 100

def read_reps() -> List[dict]:
    with REPS_TSV.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def read_members() -> Dict[str, List[str]]:
    by: Dict[str, List[str]] = defaultdict(list)
    with MEMB_TSV.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            by[r["cluster_id"]].append(r["member"])
    return by

def load_fasta(p: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    with p.open() as fh:
        acc=None; buf=[]
        for ln in fh:
            if ln.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                acc = ln[1:].strip().split()[0]
                buf=[]
            else:
                if acc is not None:
                    buf.append(ln.strip())
        if acc is not None:
            seqs[acc] = "".join(buf)
    return seqs

def shannon_entropy(seq: str) -> float:
    if not seq: return 0.0
    c = Counter(seq)
    n = len(seq)
    H = 0.0
    for k, v in c.items():
        p = v / n
        H -= p * math.log(p, 2)
    return H  # bits per symbol (max ~log2(20) ≈ 4.32)

def max_aa_fraction(seq: str) -> float:
    if not seq: return 0.0
    c = Counter(seq)
    return max(c.values()) / len(seq)

def main() -> None:
    assert REPS_TSV.exists(), f"Missing {REPS_TSV}"
    assert MEMB_TSV.exists(), f"Missing {MEMB_TSV}"
    assert FASTA26K.exists(), f"Missing {FASTA26K}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    reps = read_reps()
    memb = read_members()
    seqs = load_fasta(FASTA26K)

    # Sizes
    sizes: List[int] = []
    large_rows: List[dict] = []
    singletons = 0

    for r in reps:
        cid = r["cluster_id"]
        size = int(r.get("n_members", "0") or 0) or len(memb.get(cid, []))
        sizes.append(size)
        if size == 1:
            singletons += 1
        if size > LARGE_THRESHOLD:
            rep = r["representative"]
            s = seqs.get(rep, "")
            ent = shannon_entropy(s)
            maxfrac = max_aa_fraction(s)
            large_rows.append({
                "cluster_id": cid,
                "representative": rep,
                "size": size,
                "rep_len": len(s),
                "rep_entropy_bits": round(ent, 3),
                "rep_max_aa_fraction": round(maxfrac, 3),
                "note": "large_cluster"
            })

    # Histogram (powers of two bins up to max; plus exact bins for small sizes 1..10)
    max_size = max(sizes) if sizes else 0
    bins = []
    # small bins 1..10
    for k in range(1, 11):
        bins.append((f"{k}", sum(1 for x in sizes if x == k)))
    # geometric bins >10
    b = 16
    while b <= max(16, max_size):
        low = b//2 + 1
        high = b
        bins.append((f"{low}-{high}", sum(1 for x in sizes if low <= x <= high)))
        b *= 2
    bins.append((f">{b//2}", sum(1 for x in sizes if x > b//2)))

    # Write histogram
    with HIST_TSV.open("w") as fh:
        fh.write("bin\tcount\n")
        for name, cnt in bins:
            fh.write(f"{name}\t{cnt}\n")

    # Write large cluster table (sorted by size desc)
    large_rows.sort(key=lambda d: d["size"], reverse=True)
    with LARGE_TSV.open("w") as fh:
        fh.write("\t".join(["cluster_id","representative","size","rep_len","rep_entropy_bits","rep_max_aa_fraction","note"]) + "\n")
        for d in large_rows[:2000]:  # cap output size
            fh.write("\t".join(str(d[k]) for k in ["cluster_id","representative","size","rep_len","rep_entropy_bits","rep_max_aa_fraction","note"]) + "\n")

    # Summary
    sizes_sorted = sorted(sizes)
    def pct(p): 
        i = max(0, min(len(sizes_sorted)-1, int(p*(len(sizes_sorted)-1))))
        return sizes_sorted[i] if sizes_sorted else 0
    summary = {
        "clusters": len(reps),
        "total_members": sum(sizes),
        "singletons": singletons,
        "singleton_fraction": round(singletons / len(reps), 4) if reps else 0.0,
        "largest_cluster_size": max_size,
        "p50_cluster_size": pct(0.50),
        "p90_cluster_size": pct(0.90),
        "p99_cluster_size": pct(0.99),
        "flagged_large_clusters": len(large_rows),
        "large_threshold": LARGE_THRESHOLD
    }
    SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {HIST_TSV}")
    print(f"Wrote {LARGE_TSV}")
    print(f"Wrote {SUMMARY}")

if __name__ == "__main__":
    main()
