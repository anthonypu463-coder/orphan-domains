#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, io, json, os, sys, tempfile
from pathlib import Path
from typing import Dict, Tuple, Optional, List

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions

# Inputs (from 6A and earlier)
IDX_DIR = Path("data/structures")                       # e.g., data/structures/orphanZ70_reps/structure_index.tsv
FASTA_113K_REPS = Path("data/outputs/cluster/orphanZ70_reps.fasta")
FASTA_26K_REPS  = Path("data/outputs/cluster26k/strict26k_reps_refined.fasta")
FASTA_113K_ALL  = Path("data/outputs/orphans/orphan_zero_struct70.fasta")

# Outputs
OUT_DIR = Path("data/outputs/structures_qc")

def load_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    acc = None; buf: List[str] = []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if acc is not None: seqs[acc] = "".join(buf)
                acc = ln[1:].strip().split()[0]; buf = []
            else:
                if acc is not None: buf.append(ln.strip())
    if acc is not None: seqs[acc] = "".join(buf)
    return seqs

def clean_seq(s: str) -> str:
    return "".join(ch for ch in s.replace("\n", "").replace(" ", "") if ch.isalpha()).upper()

def mmcif_poly_seq(cif_path: Path) -> Optional[str]:
    try:
        # Biopython's MMCIF2Dict expects uncompressed; we handle .gz transparently
        if cif_path.suffix == ".gz":
            with gzip.open(cif_path, "rt") as fh:
                text = fh.read()
            d = MMCIF2Dict(io.StringIO(text))
        else:
            d = MMCIF2Dict(str(cif_path))
        for key in ("_entity_poly.pdbx_seq_one_letter_code_can", "_entity_poly.pdbx_seq_one_letter_code"):
            if key in d:
                v = d[key]
                if isinstance(v, list):
                    candidates = [clean_seq(x) for x in v if isinstance(x, str)]
                    if candidates:
                        return max(candidates, key=len)
                elif isinstance(v, str):
                    return clean_seq(v)
    except Exception:
        return None
    return None

def mmcif_atom_seq(cif_path: Path) -> Optional[str]:
    try:
        # MMCIFParser requires a filename; write a temp file if gz
        if cif_path.suffix == ".gz":
            with gzip.open(cif_path, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
                tmp.write(fin.read())
                tmp_path = Path(tmp.name)
        else:
            tmp_path = cif_path
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(tmp_path))
        model = next(structure.get_models())
        ppb = PPBuilder()
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(model))
        return seq if seq else None
    except (PDBExceptions.PDBConstructionException, Exception):
        return None
    finally:
        if cif_path.suffix == ".gz":
            try: os.unlink(tmp_path)  # type: ignore[name-defined]
            except Exception: pass

def compare_simple(a: str, b: str) -> Tuple[str, int]:
    """Return ('match'|'minor'|'major', abs_len_diff). Minor if tails differ by <=5 residues total."""
    if a == b:
        return "match", 0
    da = abs(len(a) - len(b))
    # Anchor-based tolerance: allow up to 5 residues unmatched around ends
    m = min(len(a), len(b))
    # common prefix
    i = 0
    while i < m and a[i] == b[i]:
        i += 1
    # common suffix
    j = 0
    while j < m - i and a[-(j+1)] == b[-(j+1)]:
        j += 1
    unmatched = m - (i + j)
    if unmatched <= 5 and da <= 5:
        return "minor", da
    return "major", da

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["reps_113k","reps_26k","all_113k"], default="reps_113k",
                    help="Which staged set to check (must have been staged in 6A).")
    ap.add_argument("--limit", type=int, default=0, help="Cap number of entries (0 = all).")
    ap.add_argument("--workers", type=int, default=8, help="Parsing concurrency (not used in this simple version).")
    args = ap.parse_args()

    set_name = {"reps_113k":"orphanZ70_reps","reps_26k":"strict26k_reps","all_113k":"orphanZ70_all"}[args.set]
    idx = IDX_DIR / set_name / "structure_index.tsv"
    if not idx.exists():
        sys.exit(f"Missing index: {idx} (run 6A staging first for --set {args.set})")

    # FASTA for this set
    fasta = {"reps_113k": FASTA_113K_REPS, "reps_26k": FASTA_26K_REPS, "all_113k": FASTA_113K_ALL}[args.set]
    if not fasta.exists():
        sys.exit(f"Missing FASTA for {args.set}: {fasta}")
    seqs = load_fasta(fasta)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_base = OUT_DIR / f"{set_name}_seq_struct_check"
    out_tsv  = out_base.with_suffix(".tsv")
    out_json = out_base.with_suffix(".json")
    out_bad  = OUT_DIR / f"{set_name}_seq_struct_anomalies.tsv"

    checked = 0
    poly_mismatch = 0
    atom_minor = 0
    atom_major = 0
    anomalies: List[Dict[str,str]] = []

    with idx.open() as fh, out_tsv.open("w") as out:
        rdr = csv.DictReader(fh, delimiter="\t")
        out.write("\t".join([
            "uniprot","len_fasta","len_poly","len_atom",
            "poly_match","atom_status","len_diff_atom","notes"
        ]) + "\n")
        for r in rdr:
            acc = r["uniprot"]
            cif = Path(r["model_cif"])
            fasta_seq = seqs.get(acc, "")
            if not fasta_seq:
                anomalies.append({"uniprot": acc, "issue": "missing_fasta"})
                continue

            s_poly = mmcif_poly_seq(cif) or ""
            s_atom = mmcif_atom_seq(cif) or ""

            # polymer: expect exact match
            poly_ok = int(s_poly == fasta_seq)

            # atom-derived: allow minor tail diffs
            atom_status, dlen = compare_simple(s_atom, fasta_seq) if s_atom else ("major", abs(len(fasta_seq)))
            if not poly_ok: poly_mismatch += 1
            if atom_status == "minor": atom_minor += 1
            if atom_status == "major": atom_major += 1

            note_parts = []
            if not poly_ok: note_parts.append("poly_mismatch")
            if atom_status != "match": note_parts.append(f"atom_{atom_status}")
            notes = ",".join(note_parts) if note_parts else "ok"

            out.write("\t".join([
                acc, str(len(fasta_seq)), str(len(s_poly) if s_poly else 0), str(len(s_atom) if s_atom else 0),
                str(poly_ok), atom_status, str(dlen), notes
            ]) + "\n")

            if not poly_ok or atom_status == "major":
                anomalies.append({
                    "uniprot": acc,
                    "len_fasta": str(len(fasta_seq)),
                    "len_poly": str(len(s_poly) if s_poly else 0),
                    "len_atom": str(len(s_atom) if s_atom else 0),
                    "atom_status": atom_status
                })

            checked += 1
            if args.limit and checked >= args.limit:
                break

    (out_bad).write_text("uniprot\tlen_fasta\tlen_poly\tlen_atom\tatom_status\n" + "\n".join(
        f"{a['uniprot']}\t{a['len_fasta']}\t{a['len_poly']}\t{a['len_atom']}\t{a['atom_status']}" for a in anomalies
    ))

    out_json.write_text(json.dumps({
        "set": args.set,
        "checked": checked,
        "poly_exact_mismatch": poly_mismatch,
        "atom_minor": atom_minor,
        "atom_major": atom_major,
        "anomalies_tsv": str(out_bad)
    }, indent=2))

    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_json}")
    print(f"Anomalies: poly_mismatch={poly_mismatch}, atom_minor={atom_minor}, atom_major={atom_major}")

if __name__ == "__main__":
    main()
