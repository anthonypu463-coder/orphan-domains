#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, json, os, statistics, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List, Dict

from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/afdb/index/swissprot_manifest.tsv")
    ap.add_argument("--acc-list", default="data/afdb/index/swissprot_accessions.txt")
    ap.add_argument("--pae-dir", default="data/afdb/pae")
    ap.add_argument("--out", default="data/afdb/index/reference_manifest.tsv")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--limit", type=int, default=0, help="Limit number of accessions (0 = all)")
    return ap.parse_args()

def load_manifest_map(manifest_tsv: Path) -> Dict[str, Path]:
    m: Dict[str, Path] = {}
    with manifest_tsv.open() as fh:
        next(fh)
        for line in fh:
            u, p = line.rstrip("\n").split("\t", 1)
            m[u] = Path(p)
    return m

def read_accs(path: Path, limit: int) -> List[str]:
    accs = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if limit and limit < len(accs):
        accs = accs[:limit]
    return accs

def parse_stats(cif_gz: Path, pae_json: Optional[Path]) -> Tuple[Optional[int], Dict[str, float], Optional[int], str]:
    try:
        with gzip.open(cif_gz, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
            tmp.write(fin.read())
            tmp_path = tmp.name
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", tmp_path)
    except (OSError, PDBExceptions.PDBConstructionException, Exception) as e:
        try: os.unlink(tmp_path)  # type: ignore[name-defined]
        except Exception: pass
        return None, {}, None, f"parse_error:{e.__class__.__name__}"
    finally:
        try: os.unlink(tmp_path)  # type: ignore[name-defined]
        except Exception: pass

    try:
        model = next(structure.get_models())
        ppb = PPBuilder()
        seq_len = sum(len(str(pp.get_sequence())) for pp in ppb.build_peptides(model))
        plddts: List[float] = []
        for ch in model:
            for res in ch:
                ca = res["CA"] if "CA" in res else None
                if ca is not None:
                    plddts.append(float(ca.get_bfactor()))
        if not plddts:
            plddts = [float(a.get_bfactor()) for a in model.get_atoms()]
        stats = {
            "mean": statistics.fmean(plddts) if plddts else float("nan"),
            "median": statistics.median(plddts) if plddts else float("nan"),
            "min": min(plddts) if plddts else float("nan"),
            "max": max(plddts) if plddts else float("nan"),
            "frac_ge70": (sum(v >= 70.0 for v in plddts) / len(plddts)) if plddts else float("nan"),
            "frac_ge90": (sum(v >= 90.0 for v in plddts) / len(plddts)) if plddts else float("nan"),
        }
    except Exception as e:
        return None, {}, None, f"sequence_or_bfactor_error:{e.__class__.__name__}"

    pae_dim: Optional[int] = None
    if pae_json and pae_json.exists():
        try:
            data = json.loads(pae_json.read_text())
            v = None
            for k in ("predicted_aligned_error", "pae", "pae_matrix"):
                if isinstance(data, dict) and k in data:
                    v = data[k]; break
            if v is None:
                v = data
            if isinstance(v, list) and v and isinstance(v[0], list):
                n = len(v)
                if all(isinstance(r, list) and len(r) == n for r in v):
                    pae_dim = n
        except Exception:
            pae_dim = None

    return seq_len, stats, pae_dim, "ok"

def worker_task(acc: str, cif_path: str, pae_dir: str) -> List[str]:
    cif = Path(cif_path)
    pae = Path(pae_dir) / f"{acc}.json"
    seq_len, stats, pae_dim, note = parse_stats(cif, pae if pae.exists() else None)
    if seq_len is None or not stats:
        return [acc, str(cif), "", "", "", "", "", "", "", str(int(pae.exists())), str(pae_dim or ""), note]
    return [
        acc, str(cif), str(seq_len),
        f"{stats['mean']:.3f}", f"{stats['median']:.3f}", f"{stats['min']:.3f}", f"{stats['max']:.3f}",
        f"{stats['frac_ge70']:.3f}", f"{stats['frac_ge90']:.3f}",
        str(int(pae.exists())), str(pae_dim or ""), note
    ]

def main() -> None:
    args = parse_args()
    acc_list = Path(args.acc_list)
    manifest_tsv = Path(args.manifest)
    pae_dir = Path(args.pae_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    accs = read_accs(acc_list, args.limit)
    m = load_manifest_map(manifest_tsv)

    done: set[str] = set()
    if out.exists():
        with out.open() as fh:
            next(fh, None)
            for line in fh:
                done.add(line.split("\t", 1)[0])

    new_file = not out.exists()
    fout = out.open("a")
    if new_file:
        fout.write("\t".join([
            "uniprot","cif_gz","residues",
            "plddt_mean","plddt_median","plddt_min","plddt_max",
            "plddt_frac_ge70","plddt_frac_ge90",
            "pae_present","pae_dim","note"
        ]) + "\n")

    work = [(acc, str(m[acc])) for acc in accs if acc not in done and acc in m]

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(worker_task, acc, cif, str(pae_dir)): acc for acc, cif in work}
        n_done = 0
        for fut in as_completed(futs):
            row = fut.result()
            fout.write("\t".join(row) + "\n")
            n_done += 1
            if n_done % 1000 == 0:
                fout.flush()

    fout.close()
    with out.open() as fh:
        total_rows = sum(1 for _ in fh) - 1
    print(f"Wrote {out} ({total_rows} rows)")
if __name__ == "__main__":
    main()
