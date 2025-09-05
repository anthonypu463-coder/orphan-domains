#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, io, json, os, subprocess, tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Inputs
ROOT = Path(".").resolve()
STRUCT_DIR = ROOT / "data/structures"   # e.g. data/structures/orphanZ70_reps/structure_index.tsv
BYACC_ROOT = ROOT / "data/afdb/by_acc"
PAE_GLOBAL = ROOT / "data/afdb/pae"
FETCH_SH   = ROOT / "scripts/fetch_afdb_pae.sh"  # Step-1C fetcher (parallel)

def read_index(set_name: str) -> Tuple[Path, List[dict]]:
    idx = STRUCT_DIR / set_name / "structure_index.tsv"
    if not idx.exists(): raise SystemExit(f"Missing index: {idx}")
    rows: List[dict] = []
    with idx.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr: rows.append(r)
    return idx, rows

def ensure_link(src: Path, dst: Path) -> None:
    if dst.exists(): return
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(src.resolve(), dst)

def mmcif_ca_bfactors(cif_path: Path) -> Tuple[int, float, float]:
    """Return (n_CA, min_b, max_b) from CA atoms; fallback to all atoms if no CA found."""
    tmp_path: Optional[Path] = None
    try:
        if cif_path.suffix == ".gz":
            with gzip.open(cif_path, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
                tmp.write(fin.read()); tmp_path = Path(tmp.name)
        else:
            tmp_path = cif_path
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(tmp_path))
        model = next(structure.get_models())
        vals = []
        for ch in model:
            for res in ch:
                if "CA" in res: vals.append(float(res["CA"].get_bfactor()))
        if not vals:
            vals = [float(a.get_bfactor()) for a in model.get_atoms()]
        return (len(vals), (min(vals) if vals else 0.0), (max(vals) if vals else 0.0))
    except (PDBExceptions.PDBConstructionException, Exception):
        return (0, 0.0, 0.0)
    finally:
        if cif_path.suffix == ".gz" and tmp_path and tmp_path != cif_path:
            try: os.unlink(tmp_path)
            except Exception: pass

def has_poly_seq(cif_path: Path) -> bool:
    try:
        if cif_path.suffix == ".gz":
            with gzip.open(cif_path, "rt") as fh: text = fh.read()
            d = MMCIF2Dict(io.StringIO(text))
        else:
            d = MMCIF2Dict(str(cif_path))
        for k in ("_entity_poly.pdbx_seq_one_letter_code_can", "_entity_poly.pdbx_seq_one_letter_code"):
            if k in d:
                v = d[k]
                if (isinstance(v, str) and v.strip()) or (isinstance(v, list) and any(isinstance(x, str) and x.strip() for x in v)):
                    return True
    except Exception:
        return False
    return False

def fetch_pae_for(accs: List[str], parallel: int) -> int:
    if not FETCH_SH.exists():
        print("No fetcher script found; skipping download")
        return 0
    # write a temp acc list and call fetcher with ACC_LIST override
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        for a in accs: tmp.write(a + "\n")
        tmp_path = tmp.name
    env = os.environ.copy()
    env["ACC_LIST"] = tmp_path
    env["PARALLEL"] = str(parallel)
    try:
        subprocess.run([str(FETCH_SH)], check=False, env=env, cwd=str(ROOT))
    finally:
        try: os.unlink(tmp_path)
        except Exception: pass
    # count successfully fetched
    return sum((PAE_GLOBAL / f"{a}.json").exists() for a in accs)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--fetch-missing", action="store_true", help="Attempt to fetch missing PAE JSONs")
    ap.add_argument("--parallel", type=int, default=8, help="Parallel fetch slots")
    args = ap.parse_args()

    # load staged index
    idx_path, rows = read_index(args.set)
    out_dir = ROOT / "data/outputs/structures_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_tsv = out_dir / f"{args.set}_confidence_report.tsv"
    summary_js = out_dir / f"{args.set}_confidence_summary.json"

    missing_pae: List[str] = []
    no_poly: List[str] = []
    bad_plddt: List[str] = []
    staged_root = idx_path.parent

    # scan
    with report_tsv.open("w") as out:
        out.write("\t".join(["uniprot","has_poly_seq","n_b","b_min","b_max","b_in_0_100","has_pae","pae_path"]) + "\n")
        for r in rows:
            acc = r["uniprot"]
            cif = Path(r["model_cif"])
            pae_path = Path(r["pae_json"]) if r.get("pae_json","") else None
            # polymer sequence flag
            poly_ok = has_poly_seq(cif)
            if not poly_ok: no_poly.append(acc)
            # b-factor stats
            n_b, bmin, bmax = mmcif_ca_bfactors(cif)
            b_ok = int(n_b > 0 and 0.0 <= bmin <= 100.0 and 0.0 <= bmax <= 100.0)
            if not b_ok: bad_plddt.append(acc)
            # pae presence
            has_pae = "1" if (pae_path and pae_path.exists()) else "0"
            if has_pae == "0":
                # check global store; if present, link it under by_acc
                global_json = PAE_GLOBAL / f"{acc}.json"
                if global_json.exists():
                    ensure_link(global_json, staged_root / "by_acc" / acc / "pae.json")
                    has_pae = "1"
                    pae_path = staged_root / "by_acc" / acc / "pae.json"
                else:
                    missing_pae.append(acc)
            out.write("\t".join([
                acc, str(int(poly_ok)), str(n_b), f"{bmin:.2f}", f"{bmax:.2f}",
                str(b_ok), has_pae, (pae_path.as_posix() if pae_path else "")
            ]) + "\n")

    fetched = 0
    if args.fetch_missing and missing_pae:
        print(f"Fetching {len(missing_pae)} missing PAE JSONs…")
        fetched = fetch_pae_for(missing_pae, args.parallel)
        # relink any newly fetched
        for acc in missing_pae:
            src = PAE_GLOBAL / f"{acc}.json"
            if src.exists():
                ensure_link(src, staged_root / "by_acc" / acc / "pae.json")

    # summarize
    summary = {
        "set": args.set,
        "n_index": len(rows),
        "no_poly_seq": len(no_poly),
        "bad_plddt_field": len(bad_plddt),
        "missing_pae_before_fetch": len(missing_pae),
        "fetched_pae": fetched,
        "report_tsv": report_tsv.as_posix()
    }
    summary_js.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {report_tsv}")
    print(f"Wrote {summary_js}")
    if missing_pae and not args.fetch_missing:
        print("Note: PAE missing for some entries; rerun with --fetch-missing to attempt downloads.")

if __name__ == "__main__":
    main()
