#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, os, tempfile, json
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions

# Inputs (from 6A)
STRUCT_DIR = Path("data/structures")  # contains <set_name>/structure_index.tsv

# Outputs
OUT_ROOT = Path("data/outputs/structures/features")

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps",
                    help="Which staged set to process (must exist under data/structures/).")
    ap.add_argument("--limit", type=int, default=0, help="Cap number of models to process (0 = all).")
    ap.add_argument("--overwrite", action="store_true", help="Recompute even if .npz exists.")
    return ap.parse_args()

def load_index(set_name: str) -> List[dict]:
    idx = STRUCT_DIR / set_name / "structure_index.tsv"
    if not idx.exists():
        raise SystemExit(f"Missing index: {idx} (run 6A staging first)")
    with idx.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def extract_ca_and_plddt(cif_path: Path) -> Tuple[np.ndarray, np.ndarray, List[str], List[int]]:
    """
    Returns:
      ca_xyz: (N,3) float32
      plddt: (N,) float32
      chains: list of chain IDs length N
      resseq: list of residue serials length N
    """
    tmp: Optional[Path] = None
    try:
        # Biopython MMCIFParser needs a plain file; handle gz.
        if cif_path.suffix == ".gz":
            with gzip.open(cif_path, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmpf:
                tmpf.write(fin.read())
                tmp = Path(tmpf.name)
            cif_file = str(tmp)
        else:
            cif_file = str(cif_path)

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_file)
        model = next(structure.get_models())

        ca_xyz: List[List[float]] = []
        plddt: List[float] = []
        chains: List[str] = []
        resseq: List[int] = []

        for ch in model:
            ch_id = ch.id
            for res in ch:
                # Try Cα first
                if "CA" in res:
                    at = res["CA"]
                    ca_xyz.append(list(at.get_coord()))
                    plddt.append(float(at.get_bfactor()))
                    chains.append(ch_id)
                    resseq.append(int(res.get_id()[1]))
                else:
                    # Fallback: average over atoms in residue
                    atoms = list(res.get_atoms())
                    if not atoms:
                        continue
                    coords = np.array([a.get_coord() for a in atoms], dtype=np.float32)
                    bfs    = np.array([a.get_bfactor() for a in atoms], dtype=np.float32)
                    ca_xyz.append(list(coords.mean(axis=0)))
                    plddt.append(float(bfs.mean()))
                    chains.append(ch_id)
                    resseq.append(int(res.get_id()[1]))

        if not ca_xyz:
            return np.zeros((0,3), dtype=np.float32), np.zeros((0,), dtype=np.float32), [], []

        return (
            np.asarray(ca_xyz, dtype=np.float32),
            np.asarray(plddt, dtype=np.float32),
            chains,
            resseq,
        )
    except (PDBExceptions.PDBConstructionException, Exception) as e:
        # Return empty on parse error; caller logs anomaly.
        return np.zeros((0,3), dtype=np.float32), np.zeros((0,), dtype=np.float32), [], []
    finally:
        if tmp is not None:
            try: os.unlink(tmp)
            except Exception: pass

def main() -> None:
    args = parse_args()
    rows = load_index(args.set)
    if args.limit and args.limit < len(rows):
        rows = rows[:args.limit]

    out_dir = OUT_ROOT / args.set / "by_acc"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = OUT_ROOT / args.set / "features_index.tsv"
    (OUT_ROOT / args.set).mkdir(parents=True, exist_ok=True)

    done = 0
    anomalies = 0

    # Prepare/append manifest header if new
    new_file = not manifest.exists()
    mfh = manifest.open("a")
    if new_file:
        mfh.write("\t".join([
            "uniprot","n_res","plddt_min","plddt_mean","plddt_max","npz_path"
        ]) + "\n")

    for r in rows:
        acc = r["uniprot"]
        cif = Path(r["model_cif"])
        npz_path = out_dir / f"{acc}.npz"
        if npz_path.exists() and not args.overwrite:
            done += 1
            continue

        ca_xyz, plddt, chains, resseq = extract_ca_and_plddt(cif)
        if ca_xyz.shape[0] == 0:
            anomalies += 1
            # Write an empty stub for reproducibility
            np.savez_compressed(npz_path, acc=acc, ca_xyz=np.zeros((0,3), np.float32),
                                plddt=np.zeros((0,), np.float32), chains=np.array([], object),
                                resseq=np.array([], np.int32))
            mfh.write(f"{acc}\t0\t0.0\t0.0\t0.0\t{npz_path.as_posix()}\n")
            done += 1
            continue

        np.savez_compressed(
            npz_path,
            acc=acc,
            ca_xyz=ca_xyz,
            plddt=plddt.astype(np.float32),
            chains=np.array(chains, dtype=object),
            resseq=np.asarray(resseq, dtype=np.int32),
        )
        mfh.write(
            f"{acc}\t{ca_xyz.shape[0]}\t{plddt.min():.2f}\t{plddt.mean():.2f}\t{plddt.max():.2f}\t{npz_path.as_posix()}\n"
        )
        done += 1

    mfh.close()
    # Write a tiny summary sidecar
    summary = {
        "set": args.set,
        "processed": done,
        "anomalies": anomalies,
        "features_dir": str(out_dir),
        "manifest": manifest.as_posix(),
    }
    (OUT_ROOT / args.set / "features_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
