#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, os, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, Dict, List

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions

def clean_seq(s: str) -> str:
    return "".join(ch for ch in s.replace("\n", "").replace(" ", "") if ch.isalpha()).upper()

def seq_from_mmcif_dict(cif_path: Path) -> Optional[str]:
    d = MMCIF2Dict(str(cif_path))
    for k in ("_entity_poly.pdbx_seq_one_letter_code_can", "_entity_poly.pdbx_seq_one_letter_code"):
        if k in d:
            v = d[k]
            if isinstance(v, list):
                candidates = [clean_seq(x) for x in v if isinstance(x, str)]
                if candidates:
                    return max(candidates, key=len)
            elif isinstance(v, str):
                return clean_seq(v)
    return None

def seq_from_atoms(cif_path: Path) -> Optional[str]:
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(cif_path))
        model = next(structure.get_models())
        ppb = PPBuilder()
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(model))
        return seq if seq else None
    except (PDBExceptions.PDBConstructionException, Exception):
        return None

def read_seq_from_cif_gz(cif_gz: Path) -> Optional[str]:
    with gzip.open(cif_gz, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
        tmp.write(fin.read())
        tmp_path = Path(tmp.name)
    try:
        seq = seq_from_mmcif_dict(tmp_path)
        if not seq:
            seq = seq_from_atoms(tmp_path)
        return seq
    finally:
        try: tmp_path.unlink()
        except Exception: pass

def load_manifest(manifest_tsv: Path) -> Dict[str, Path]:
    m: Dict[str, Path] = {}
    with manifest_tsv.open() as fh:
        header = next(fh, None)
        for line in fh:
            acc, p = line.rstrip("\n").split("\t", 1)
            m[acc] = Path(p)
    return m

def task(acc: str, cif_gz: Path) -> Tuple[str, Optional[str]]:
    try:
        return acc, read_seq_from_cif_gz(cif_gz)
    except Exception:
        return acc, None

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/afdb/index/swissprot_manifest.tsv")
    ap.add_argument("--acc-list", default="data/afdb/index/swissprot_accessions.txt")
    ap.add_argument("--out-fasta", default="data/afdb/index/swissprot_sequences.fasta")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    man = load_manifest(Path(args.manifest))
    accs = [ln.strip() for ln in Path(args.acc_list).read_text().splitlines() if ln.strip()]
    if args.limit and args.limit < len(accs):
        accs = accs[:args.limit]

    out = Path(args.out_fasta)
    out.parent.mkdir(parents=True, exist_ok=True)

    done: set[str] = set()
    if out.exists():
        with out.open() as fh:
            for ln in fh:
                if ln.startswith(">"):
                    done.add(ln[1:].strip())
    todo = [(a, man[a]) for a in accs if a not in done and a in man]

    # Append mode to allow resume
    fout = out.open("a")
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(task, a, p): a for a, p in todo}
        for fut in as_completed(futs):
            acc, seq = fut.result()
            if not seq:
                continue
            fout.write(f">{acc}\n")
            # wrap to 60 cols
            for i in range(0, len(seq), 60):
                fout.write(seq[i:i+60] + "\n")
    fout.close()
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
