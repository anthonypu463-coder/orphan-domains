#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, os
from pathlib import Path
from typing import Dict, List, Tuple

# Inputs from earlier steps
MANIFEST   = Path("data/afdb/index/reference_manifest.tsv")         # uniprot, residues, plddt_mean
BYACC_ROOT = Path("data/afdb/by_acc")                               # per-ACC symlinks (model_v4.cif.gz, pae.json)
CIF_MAN    = Path("data/afdb/index/swissprot_manifest.tsv")         # fallback mmCIF path lookup
PAE_DIR    = Path("data/afdb/pae")                                  # fallback PAE json

# Representative sets / universes
REPS_113K  = Path("data/outputs/cluster/orphanZ70_reps.fasta")      # default (publication primary, non-redundant)
REPS_26K   = Path("data/outputs/cluster26k/strict26k_reps_refined.fasta")
ALL_113K   = Path("data/outputs/orphans/orphan_zero_struct70_accessions.txt")

def read_ids_from_fasta(fa: Path) -> List[str]:
    ids: List[str] = []
    with fa.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                ids.append(ln[1:].strip().split()[0])
    return ids

def read_ids_from_list(txt: Path) -> List[str]:
    return [ln.strip() for ln in txt.read_text().splitlines() if ln.strip()]

def load_manifest(path: Path) -> Tuple[Dict[str,int], Dict[str,float]]:
    L: Dict[str,int] = {}
    P: Dict[str,float] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"].strip()
            if not acc: continue
            try:
                L[acc] = int(r.get("residues","") or 0)
            except: L[acc] = 0
            try:
                P[acc] = float(r.get("plddt_mean","") or 0.0)
            except: P[acc] = 0.0
    return L, P

def load_cif_map(path: Path) -> Dict[str,Path]:
    m: Dict[str,Path] = {}
    with path.open() as fh:
        next(fh, None)  # header
        for ln in fh:
            acc, cif = ln.rstrip("\n").split("\t", 1)
            m[acc] = Path(cif)
    return m

def ensure_symlink(src: Path, dst: Path) -> None:
    if dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    # Use absolute source paths for robustness
    os.symlink(src.resolve(), dst)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["reps_113k","reps_26k","all_113k"], default="reps_113k",
                    help="Which proteins to stage for Step 6 analysis.")
    ap.add_argument("--limit", type=int, default=0, help="Optional cap for quick smoke tests.")
    args = ap.parse_args()

    # Pick ID source and workspace name
    if args.set == "reps_113k":
        ids = read_ids_from_fasta(REPS_113K); set_name = "orphanZ70_reps"
    elif args.set == "reps_26k":
        ids = read_ids_from_fasta(REPS_26K);  set_name = "strict26k_reps"
    else:
        ids = read_ids_from_list(ALL_113K);   set_name = "orphanZ70_all"

    if args.limit and args.limit < len(ids):
        ids = ids[:args.limit]

    out_root = Path("data/structures") / set_name
    index_tsv = out_root / "structure_index.tsv"
    out_root.mkdir(parents=True, exist_ok=True)

    # Support lookup
    lengths, plddts = load_manifest(MANIFEST)
    cif_map = load_cif_map(CIF_MAN)

    staged = 0
    missing_cif = 0
    pae_missing = 0

    with index_tsv.open("w") as out:
        out.write("\t".join([
            "uniprot","residues","plddt_mean",
            "model_cif","pae_json","has_pae","model_bytes"
        ]) + "\n")

        for acc in ids:
            # Prefer by_acc layout from Step 1E
            by_acc_dir = BYACC_ROOT / acc
            model = by_acc_dir / "model_v4.cif.gz"
            pae   = by_acc_dir / "pae.json"

            if not model.exists():
                # Fallback to manifest path
                cif_src = cif_map.get(acc)
                if cif_src and Path(cif_src).exists():
                    model = cif_src
                else:
                    missing_cif += 1
                    continue

            # If by_acc dir exists but symlinks are missing, create them for a clean workspace
            dest_dir = out_root / "by_acc" / acc
            ensure_symlink(Path(model), dest_dir / "model_v4.cif.gz")

            # PAE resolution and link (best-effort)
            pae_src = nae = None
            if (BYACC_ROOT / acc / "pae.json").exists():
                pae_src = BYACC_ROOT / acc / "pae.json"
            elif (PAE_DIR / f"{acc}.json").exists():
                pae_src = PAE_DIR / f"{acc}.json"
            if pae_src:
                ensure_symlink(pae_src, dest_dir / "pae.json")
                has_pae = "1"
                pae_path = (dest_dir / "pae.json").as_posix()
            else:
                has_pae = "0"
                pae_path = ""

            # Stats
            L = lengths.get(acc, 0)
            P = plddts.get(acc, 0.0)
            model_bytes = (dest_dir / "model_v4.cif.gz").stat().st_size if (dest_dir / "model_v4.cif.gz").exists() else 0

            out.write("\t".join([
                acc, str(L), f"{P:.2f}",
                (dest_dir / "model_v4.cif.gz").as_posix(),
                pae_path, has_pae, str(model_bytes)
            ]) + "\n")
            staged += 1
            if has_pae == "0":
                pae_missing += 1

    print(f"Workspace: {out_root}")
    print(f"Index:     {index_tsv}")
    print(f"Staged models: {staged} | Missing mmCIF: {missing_cif} | Missing PAE: {pae_missing}")

if __name__ == "__main__":
    main()
