#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
from collections import defaultdict, Counter

REPS_TSV = Path("data/outputs/cluster26k/cluster26k_representatives.tsv")  # cluster_id, representative, n_members, ...
MEMB_TSV = Path("data/outputs/cluster26k/cluster26k_members.tsv")          # cluster_id, member

ASSIGN_TSV = Path("data/outputs/cluster26k/strict26k_assignments.tsv")     # uniprot, cluster_id, representative
SIZE_TSV   = Path("data/outputs/cluster26k/strict26k_size_counts.tsv")     # size,count
SUMMARY_JS = Path("data/outputs/cluster26k/strict26k_cluster_summary.json")

def main() -> None:
    assert REPS_TSV.exists(), f"Missing {REPS_TSV}"
    assert MEMB_TSV.exists(), f"Missing {MEMB_TSV}"

    # Map cluster → representative
    cid_to_rep: dict[str, str] = {}
    with REPS_TSV.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            cid_to_rep[r["cluster_id"]] = r["representative"]

    # Collect members
    cid_to_members: dict[str, list[str]] = defaultdict(list)
    with MEMB_TSV.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            cid_to_members[r["cluster_id"]].append(r["member"])

    # Assignment file
    with ASSIGN_TSV.open("w") as fh:
        fh.write("\t".join(["uniprot","cluster_id","representative"]) + "\n")
        for cid, mems in cid_to_members.items():
            rep = cid_to_rep.get(cid, "")
            for m in mems:
                fh.write(f"{m}\t{cid}\t{rep}\n")

    # Size distribution
    sizes = [len(m) for m in cid_to_members.values()]
    size_counts = Counter(sizes)
    with SIZE_TSV.open("w") as fh:
        fh.write("size\tcount\n")
        for s in sorted(size_counts):
            fh.write(f"{s}\t{size_counts[s]}\n")

    # Summary
    clusters = len(cid_to_members)
    members  = sum(sizes)
    singletons = size_counts.get(1, 0)
    summary = {
        "clusters": clusters,
        "total_members": members,
        "singletons": singletons,
        "singleton_fraction": round(singletons / clusters, 6) if clusters else 0.0,
        "largest_cluster_size": max(sizes) if sizes else 0,
        "size_counts": dict(sorted(size_counts.items()))
    }
    SUMMARY_JS.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {ASSIGN_TSV}")
    print(f"Wrote {SIZE_TSV}")
    print(f"Wrote {SUMMARY_JS}")

if __name__ == "__main__":
    main()
