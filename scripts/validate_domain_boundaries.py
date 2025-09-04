#!/usr/bin/env python3
from __future__ import annotations
import csv, json, random
from pathlib import Path
from typing import Dict, List, Tuple, Optional

FINAL = Path("data/outputs/architecture/pfam_domains_final.tsv")          # 3A
CONF  = Path("data/outputs/hmmscan/pfam_hits_confirmed.tsv")              # 2B
DROP2C= Path("data/outputs/hmmscan/pfam_domains_dropped.tsv")             # 2C
OUTDIR= Path("data/outputs/architecture/validation")
SAMPLE_TSV = OUTDIR / "domain_boundary_sample.tsv"
OVERLAP_TSV= OUTDIR / "overlap_decision_checks.tsv"
SUMMARY_JSON= OUTDIR / "summary.json"

def to_int(s: str) -> int:
    try: return int(s)
    except Exception: return 0
def to_float(s: str) -> Optional[float]:
    try: return float(s)
    except Exception: return None
def norm_acc(acc: str) -> str:
    acc = (acc or "").strip()
    i = acc.find(".")
    return acc[:i] if i > 0 else acc

def load_final(path: Path) -> List[dict]:
    rows: List[dict] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            rows.append({
                "uniprot": r["uniprot"],
                "pfam_acc": r["pfam_acc"],
                "pfam_id":  r["pfam_id"],
                "start":    to_int(r["start"]),
                "end":      to_int(r["end"]),
                "dom_score":to_float(r.get("dom_score","")),
                "full_score":to_float(r.get("full_score","")),
                "note":     r.get("note",""),
            })
    return rows

def load_confirmed(path: Path) -> Dict[Tuple[str,str], List[dict]]:
    by: Dict[Tuple[str,str], List[dict]] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            k = (r["uniprot"], norm_acc(r["pfam_acc"] or r["pfam_id"]))
            by.setdefault(k, []).append({
                "env_from": to_int(r["env_from"]), "env_to": to_int(r["env_to"]),
                "hmm_from": to_int(r["hmm_from"]), "hmm_to": to_int(r["hmm_to"]),
                "tlen":     to_int(r["tlen"]),
                "dom_score":to_float(r.get("dom_score","")), "full_score": to_float(r.get("full_score","")),
            })
    # sort each list by dom_score desc, then longer env
    for lst in by.values():
        lst.sort(key=lambda x: (-(x["dom_score"] or float("-inf")), -(x["env_to"]-x["env_from"]+1)))
    return by

def load_dropped(path: Path) -> List[dict]:
    if not path.exists(): return []
    rows: List[dict] = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            rows.append({
                "uniprot": r["uniprot"],
                "pfam_acc": r["pfam_acc"],
                "pfam_id":  r["pfam_id"],
                "env_from": to_int(r["env_from"]),
                "env_to":   to_int(r["env_to"]),
                "dom_score":to_float(r.get("dom_score","")),
                "full_score":to_float(r.get("full_score","")),
                "reason":   r.get("reason",""),
            })
    return rows

def ovlen(a1,a2,b1,b2)->int:
    lo,hi = max(a1,b1), min(a2,b2)
    return max(0, hi-lo+1)

def best_matching_hit(conf_hits: List[dict], s:int, e:int) -> Optional[dict]:
    best=None; best_ov=-1
    for h in conf_hits:
        o = ovlen(s,e,h["env_from"],h["env_to"])
        if o > best_ov:
            best, best_ov = h, o
    return best

def main(sample_n: int = 200, seed: int = 7):
    if not FINAL.exists() or not CONF.exists():
        raise SystemExit("Missing inputs; ensure 3A and 2B completed.")
    OUTDIR.mkdir(parents=True, exist_ok=True)

    final_rows = load_final(FINAL)
    confirmed  = load_confirmed(CONF)
    dropped    = load_dropped(DROP2C)

    # Sample domain rows
    random.seed(seed)
    sample = final_rows if len(final_rows) <= sample_n else random.sample(final_rows, sample_n)

    short_flags = 0
    trimmed_flags=0
    with SAMPLE_TSV.open("w") as fh:
        fh.write("\t".join([
            "uniprot","pfam_acc","pfam_id","start","end","final_len",
            "tlen","frac_len_vs_hmm","orig_env_from","orig_env_to","orig_env_len",
            "trimmed","dom_score_final","full_score_final","hmm_cov_orig","notes"
        ]) + "\n")
        for r in sample:
            key = (r["uniprot"], norm_acc(r["pfam_acc"] or r["pfam_id"]))
            hits = confirmed.get(key, [])
            m = best_matching_hit(hits, r["start"], r["end"]) if hits else None
            final_len = r["end"] - r["start"] + 1
            tlen = (m["tlen"] if m else 0)
            frac_vs_hmm = (final_len / tlen) if tlen > 0 else 0.0
            orig_env_len = (m["env_to"] - m["env_from"] + 1) if m else 0
            trimmed = int(m is not None and final_len < orig_env_len)
            hmm_cov = ((m["hmm_to"] - m["hmm_from"] + 1) / tlen) if (m and tlen>0) else 0.0
            notes = []
            if tlen>0 and (final_len < 0.4*tlen or final_len < 30):
                notes.append("short_vs_hmm")
                short_flags += 1
            if trimmed:
                trimmed_flags += 1
            fh.write("\t".join([
                r["uniprot"], r["pfam_acc"], r["pfam_id"], str(r["start"]), str(r["end"]), str(final_len),
                str(tlen), f"{frac_vs_hmm:.3f}", str(m["env_from"] if m else ""), str(m["env_to"] if m else ""), str(orig_env_len),
                str(trimmed), f"{(r['dom_score'] or 0):.2f}", f"{(r['full_score'] or 0):.2f}", f"{hmm_cov:.3f}",
                ",".join(notes) if notes else "ok"
            ]) + "\n")

    # Overlap decision checks: dropped vs kept score ordering
    kept_by_acc: Dict[str, List[Tuple[int,int,Optional[float],Optional[float],str]]] = {}
    for r in final_rows:
        kept_by_acc.setdefault(r["uniprot"], []).append((r["start"], r["end"], r["dom_score"], r["full_score"], r["pfam_acc"]))
    weaker_ok=0; weaker_bad=0
    with OVERLAP_TSV.open("w") as fh:
        fh.write("\t".join(["uniprot","dropped_acc","kept_acc","overlap","drop_dom_score","keep_dom_score","ok"]) + "\n")
        for d in dropped:
            if d.get("reason","") not in ("overlap_drop",""):  # focus on overlap drops
                continue
            kept_iv = kept_by_acc.get(d["uniprot"], [])
            best=None; best_o=0
            for s,e,ds,fs,kacc in kept_iv:
                o = ovlen(d["env_from"], d["env_to"], s, e)
                if o>best_o:
                    best=(s,e,ds,fs,kacc); best_o=o
            if not best or best_o==0:
                continue
            _,_,k_dom,_,kacc = best
            kd = k_dom if k_dom is not None else float("-inf")
            dd = d["dom_score"] if d["dom_score"] is not None else float("-inf")
            ok = int(kd >= dd)
            if ok: weaker_ok += 1
            else: weaker_bad += 1
            fh.write("\t".join([
                d["uniprot"], d["pfam_acc"], kacc, str(best_o),
                f"{(d['dom_score'] or 0):.2f}", f"{(k_dom or 0):.2f}", str(ok)
            ]) + "\n")

    SUMMARY_JSON.write_text(json.dumps({
        "sample_size": len(sample),
        "short_vs_hmm_flags": short_flags,
        "trimmed_flags": trimmed_flags,
        "overlap_checks_ok": weaker_ok,
        "overlap_checks_bad": weaker_bad
    }, indent=2))
    print(f"Wrote {SAMPLE_TSV}")
    print(f"Wrote {OVERLAP_TSV}")
    print(f"Wrote {SUMMARY_JSON}")

if __name__ == "__main__":
    main()
