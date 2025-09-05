#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, gzip, os, sys, tempfile, shutil
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.PDB import MMCIFParser, PDBIO, Select, PDBExceptions

SEG_ROOT   = Path("data/outputs/segments")                 # from 7D: segments_final.tsv
STRUCT_DIR = Path("data/structures")                       # from 6A: by_acc/<ACC>/model_v4.cif.gz

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["orphanZ70_reps","strict26k_reps","orphanZ70_all"], default="orphanZ70_reps")
    ap.add_argument("--limit", type=int, default=0, help="limit proteins (0=all)")
    ap.add_argument("--compress", action="store_true", default=True, help="write .pdb.gz instead of .pdb")
    return ap.parse_args()

def load_segments(set_name: str) -> Dict[str, List[Tuple[str,int,int,int]]]:
    seg_f = SEG_ROOT / set_name / "segments_final.tsv"
    if not seg_f.exists():
        raise SystemExit(f"Missing {seg_f}. Run 7D first.")
    by_acc: Dict[str, List[Tuple[str,int,int,int]]] = {}
    with seg_f.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            acc = r["uniprot"]
            fid = r["fragment_id"]
            s,e,L = int(r["start"]), int(r["end"]), int(r["length"])
            by_acc.setdefault(acc, []).append((fid, s, e, L))
    for acc in by_acc:
        by_acc[acc].sort(key=lambda x: (x[1], x[2]))
    return by_acc

def staged_cif_path(set_name: str, acc: str) -> Path:
    return STRUCT_DIR / set_name / "by_acc" / acc / "model_v4.cif.gz"

class ResidueRangeSelector(Select):
    def __init__(self, start: int, end: int):
        super().__init__()
        self.start = start
        self.end = end
    def accept_residue(self, residue):
        het, resseq, icode = residue.get_id()
        if het.strip():    # skip hetero/non-standard residues
            return 0
        return 1 if (self.start <= int(resseq) <= self.end) else 0

def write_segment_pdb(cif_gz: Path, out_pdb_gz: Path, start: int, end: int, compress: bool) -> bool:
    tmp_cif = None
    tmp_pdb = None
    try:
        # Decompress CIF.GZ to a temp CIF file for parsing
        with gzip.open(cif_gz, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tf:
            shutil.copyfileobj(fin, tf)
            tmp_cif = Path(tf.name)

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(tmp_cif))

        # Save slice to a temp PDB (plain text), then gzip if requested
        io = PDBIO()
        io.set_structure(structure)
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as pf:
            tmp_pdb = Path(pf.name)
        io.save(str(tmp_pdb), select=ResidueRangeSelector(start, end))

        if compress:
            with open(tmp_pdb, "rb") as fin, gzip.open(out_pdb_gz, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            # Write plain .pdb (rename target path sans .gz)
            out_pdb = out_pdb_gz.with_suffix("") if out_pdb_gz.suffix == ".gz" else out_pdb_gz
            shutil.move(str(tmp_pdb), str(out_pdb))
            tmp_pdb = None  # moved

        return out_pdb_gz.exists() and out_pdb_gz.stat().st_size > 0
    except (PDBExceptions.PDBConstructionException, Exception):
        return False
    finally:
        for p in (tmp_cif, tmp_pdb):
            try:
                if p and Path(p).exists():
                    Path(p).unlink()
            except Exception:
                pass

def main() -> None:
    args = parse_args()
    set_name = args.set

    segments = load_segments(set_name)
    if args.limit and args.limit < len(segments):
        keep = dict(sorted(segments.items())[:args.limit])
        segments = keep

    out_dir = SEG_ROOT / set_name
    pdb_dir = out_dir / "pdb"
    pdb_dir.mkdir(parents=True, exist_ok=True)

    tsv = out_dir / "segments_table.tsv"
    jsl = out_dir / "segments.jsonl"
    summary = out_dir / "segments_output_summary.json"

    total_frags = ok_frags = fail_frags = prots = 0

    with tsv.open("w") as th, jsl.open("w") as jh:
        th.write("\t".join(["uniprot","fragment_id","start","end","length","pdb_path"]) + "\n")

        for acc, frags in segments.items():
            prots += 1
            cif = staged_cif_path(set_name, acc)
            if not cif.exists():
                continue
            for (fid, s, e, L) in frags:
                total_frags += 1
                out_name = f"{acc}_{fid}.pdb.gz" if args.compress else f"{acc}_{fid}.pdb"
                out_pdb = pdb_dir / out_name
                if write_segment_pdb(cif, out_pdb, s, e, args.compress):
                    ok_frags += 1
                    th.write("\t".join([acc, fid, str(s), str(e), str(L), out_pdb.as_posix()]) + "\n")
                    jh.write(json.dumps({
                        "uniprot": acc, "fragment_id": fid, "start": s, "end": e, "length": L,
                        "pdb_path": out_pdb.as_posix()
                    }) + "\n")
                else:
                    fail_frags += 1
                    try:
                        if out_pdb.exists(): out_pdb.unlink()
                    except Exception:
                        pass

    summary.write_text(json.dumps({
        "set": set_name,
        "proteins": prots,
        "fragments_total": total_frags,
        "fragments_ok": ok_frags,
        "fragments_failed": fail_frags,
        "pdb_dir": pdb_dir.as_posix(),
        "segments_table": tsv.as_posix(),
        "segments_jsonl": jsl.as_posix()
    }, indent=2))
    print(f"Wrote {tsv}")
    print(f"Wrote {jsl}")
    print(f"Wrote {summary}")

if __name__ == "__main__":
    main()
