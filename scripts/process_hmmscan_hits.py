#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
from typing import Dict, Tuple, List, Optional

def as_int(s: str) -> Optional[int]:
    try: return int(s)
    except Exception: return None

def as_float(s: str) -> Optional[float]:
    try: return float(s)
    except Exception: return None

def load_ga_table(path: Path) -> Tuple[Dict[str, Tuple[float,float]], Dict[str, Tuple[float,float]]]:
    by_acc: Dict[str, Tuple[float,float]] = {}
    by_id:  Dict[str, Tuple[float,float]] = {}
    with path.open() as fh:
        next(fh)
        for ln in fh:
            pfam_acc, pfam_id, ga_seq, ga_dom = ln.rstrip("\n").split("\t")
            by_acc[pfam_acc] = (float(ga_seq), float(ga_dom))
            by_id[pfam_id]   = (float(ga_seq), float(ga_dom))
    return by_acc, by_id

def parse_domtbl_line(ln: str):
    # domtblout fields (hmmscan): see HMMER3 manual.
    p = ln.strip().split()
    if len(p) < 23:  # defensive
        return None
    tname, tacc = p[0], (p[1] if p[1] != "-" else "")
    tlen = as_int(p[2])
    qname, qlen = p[3], as_int(p[5])
    full_e = p[6]
    full_s = as_float(p[7])
    i_eval = p[12]
    dscore = as_float(p[13])
    hmm_from, hmm_to = as_int(p[15]), as_int(p[16])
    env_from, env_to = as_int(p[19]), as_int(p[20])
    desc = " ".join(p[22:]) if len(p) > 22 else ""
    dom_len = (env_to - env_from + 1) if (env_from is not None and env_to is not None) else None
    seq_cov = (dom_len / qlen) if (dom_len is not None and qlen and qlen > 0) else None
    hmm_cov = ((hmm_to - hmm_from + 1) / tlen) if (hmm_from is not None and hmm_to is not None and tlen and tlen > 0) else None
    return {
        "uniprot": qname,
        "pfam_id": tname,
        "pfam_acc": tacc,
        "qlen": qlen or 0,
        "tlen": tlen or 0,
        "full_evalue": full_e,
        "full_score": full_s,
        "i_evalue": i_eval,
        "dom_score": dscore,
        "hmm_from": hmm_from or 0,
        "hmm_to": hmm_to or 0,
        "env_from": env_from or 0,
        "env_to": env_to or 0,
        "dom_len": dom_len or 0,
        "seq_cov_frac": f"{seq_cov:.4f}" if seq_cov is not None else "",
        "hmm_cov_frac": f"{hmm_cov:.4f}" if hmm_cov is not None else "",
        "description": desc,
    }

def classify(rec, ga_by_acc, ga_by_id):
    ga = None
    if rec["pfam_acc"] and rec["pfam_acc"] in ga_by_acc:
        ga = ga_by_acc[rec["pfam_acc"]]
    elif rec["pfam_id"] in ga_by_id:
        ga = ga_by_id[rec["pfam_id"]]
    if ga is None:
        return "unknown"
    ga_seq, ga_dom = ga
    fs, ds = rec["full_score"], rec["dom_score"]
    if fs is None and ds is None:
        return "unknown"
    if (fs is not None and fs >= ga_seq) or (ds is not None and ds >= ga_dom):
        return "confirmed"
    return "sub_ga"

def main() -> None:
    root = Path("data/outputs/hmmscan")
    src = root / "pfam_domtblout.all"
    ga_path = root / "pfam_ga.tsv"
    if not src.exists(): sys.exit(f"Missing {src}")
    if not ga_path.exists(): sys.exit(f"Missing {ga_path} (run build_pfam_ga_table.py)")

    ga_by_acc, ga_by_id = load_ga_table(ga_path)

    records: List[dict] = []
    with src.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            rec = parse_domtbl_line(ln)
            if rec: records.append(rec)

    for r in records:
        r["class"] = classify(r, ga_by_acc, ga_by_id)

    records.sort(key=lambda r: (r["uniprot"], r["env_from"], r["env_to"]))

    cols = ["uniprot","pfam_id","pfam_acc","qlen","tlen","full_evalue","full_score","i_evalue","dom_score",
            "hmm_from","hmm_to","env_from","env_to","dom_len","seq_cov_frac","hmm_cov_frac","class","description"]
    out_conf = root / "pfam_hits_confirmed.tsv"
    out_sub  = root / "pfam_hits_subga.tsv"
    out_unk  = root / "pfam_hits_unknown.tsv"

    def write_filtered(path: Path, keep: set[str]) -> int:
        n = 0
        with path.open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in records:
                if r["class"] in keep:
                    fh.write("\t".join(str(r[c]) for c in cols) + "\n")
                    n += 1
        return n

    n_conf = write_filtered(out_conf, {"confirmed"})
    n_sub  = write_filtered(out_sub, {"sub_ga"})
    n_unk  = write_filtered(out_unk, {"unknown"})
    print(f"Wrote {out_conf} ({n_conf} rows)")
    print(f"Wrote {out_sub} ({n_sub} rows)")
    print(f"Wrote {out_unk} ({n_unk} rows)")

if __name__ == "__main__":
    main()
