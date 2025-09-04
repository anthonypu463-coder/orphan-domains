#!/usr/bin/env python3
from pathlib import Path

IN = Path("data/outputs/hmmscan/pfam_domtblout.all")
OUT = Path("data/outputs/hmmscan/pfam_hits.tsv")

def is_int(s: str) -> bool:
    try:
        int(s); return True
    except: return False

def parse_line(line: str):
    s = line.strip()
    if not s or s.startswith("#"):
        return None
    p = s.split()
    # Require the core 22+ columns; guard numeric fields to avoid stray text lines
    if len(p) < 22:  # description may be absent; valid data lines have >=22 tokens
        return None
    if not (is_int(p[2]) and is_int(p[5]) and is_int(p[15]) and is_int(p[16]) and is_int(p[19]) and is_int(p[20])):
        return None

    t_name = p[0]
    t_acc  = "" if p[1] == "-" else p[1]
    tlen   = int(p[2])
    q_name = p[3]
    qlen   = int(p[5])
    i_eval = p[12]
    dscore = p[13]
    hmm_from, hmm_to = int(p[15]), int(p[16])
    env_from, env_to = int(p[19]), int(p[20])
    desc = " ".join(p[22:]) if len(p) > 22 else ""

    dom_len = env_to - env_from + 1
    seq_cov = f"{(dom_len/qlen):.4f}" if qlen > 0 else ""
    hmm_cov = f"{((hmm_to-hmm_from+1)/tlen):.4f}" if tlen > 0 else ""

    return [q_name, t_name, t_acc, str(tlen), str(qlen), i_eval, dscore,
            str(hmm_from), str(hmm_to), str(env_from), str(env_to),
            str(dom_len), seq_cov, hmm_cov, desc]

def main():
    if not IN.exists():
        raise SystemExit(f"Missing input: {IN}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    kept = skipped = 0
    with IN.open() as fin, OUT.open("w") as fout:
        fout.write("\t".join([
            "uniprot","pfam_name","pfam_acc","tlen","qlen","i_evalue","score",
            "hmm_from","hmm_to","env_from","env_to","dom_len","seq_cov_frac","hmm_cov_frac","description"
        ]) + "\n")
        for ln in fin:
            rec = parse_line(ln)
            if rec:
                fout.write("\t".join(rec) + "\n"); kept += 1
            else:
                skipped += 1
    print(f"Wrote {OUT} (kept={kept}, skipped={skipped})")

if __name__ == "__main__":
    main()
