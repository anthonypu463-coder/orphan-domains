#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, io, json, os, random, statistics, sys, tempfile, time, urllib.request, urllib.error
from pathlib import Path
from typing import Optional, Tuple, List, Any

try:
    from orphan.logging_utils import configure_logging
    import logging
    configure_logging(level="INFO", json_format=True, quiet=False)
except Exception:
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
log = logging.getLogger("verify.integrity")

from Bio.PDB import MMCIFParser, PPBuilder, PDBExceptions

def read_accessions(acc_path: Path) -> List[str]:
    accs = [line.strip() for line in acc_path.read_text().splitlines() if line.strip()]
    if not accs:
        raise SystemExit(f"No accessions found in {acc_path}")
    return accs

def find_cif_gz_for(acc: str, manifest_tsv: Path) -> Optional[Path]:
    with manifest_tsv.open() as fh:
        next(fh)
        for line in fh:
            u, p = line.rstrip("\n").split("\t", 1)
            if u == acc:
                return Path(p)
    return None

def parse_structure_from_gz(cif_gz: Path) -> Tuple[Optional[object], Optional[int], List[float], str]:
    try:
        with gzip.open(cif_gz, "rb") as fin, tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
            tmp.write(fin.read())
            tmp_path = tmp.name
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", tmp_path)
    except (OSError, PDBExceptions.PDBConstructionException, Exception) as e:
        return None, None, [], f"parse_error:{e.__class__.__name__}"
    finally:
        try: os.unlink(tmp_path)  # type: ignore[name-defined]
        except Exception: pass

    ppb = PPBuilder()
    try:
        model = next(structure.get_models())
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(model))
        b_factors = [float(at.get_bfactor()) for at in model.get_atoms()]
        return structure, len(seq), b_factors, "ok"
    except Exception as e:
        return structure, None, [], f"sequence_build_error:{e.__class__.__name__}"

def _find_square_numeric_matrix(obj: Any, max_depth: int = 5) -> Optional[int]:
    """Return side length if a square matrix of numbers is found anywhere (nested)."""
    if max_depth < 0:
        return None
    # direct list-of-lists
    if isinstance(obj, list) and obj and isinstance(obj[0], list):
        n = len(obj)
        for row in obj:
            if not isinstance(row, list) or len(row) != n:
                return None
            for x in row:
                if not isinstance(x, (int, float)):
                    return None
        return n
    # dict: search values
    if isinstance(obj, dict):
        for v in obj.values():
            n = _find_square_numeric_matrix(v, max_depth - 1)
            if n is not None:
                return n
    # list: search elements
    if isinstance(obj, list):
        for v in obj:
            n = _find_square_numeric_matrix(v, max_depth - 1)
            if n is not None:
                return n
    return None

def load_pae_json(pae_json: Path) -> Tuple[bool, Optional[int], str]:
    if not pae_json.exists():
        return False, None, "missing_pae"
    try:
        data = json.loads(pae_json.read_text())
    except Exception as e:
        return True, None, f"json_load_error:{e.__class__.__name__}"
    n = None
    # Try common keys first for speed
    for k in ("predicted_aligned_error", "pae", "pae_matrix"):
        v = data.get(k) if isinstance(data, dict) else None
        if v is not None:
            n = _find_square_numeric_matrix(v, max_depth=2)
            if n is not None:
                return True, n, "ok"
    # Fallback: search entire structure
    n = _find_square_numeric_matrix(data, max_depth=5)
    if n is not None:
        return True, n, "ok"
    return True, None, "no_square_matrix_found"

def fetch_uniprot_fasta(acc: str, timeout: int = 20, retries: int = 3) -> Optional[str]:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    req = urllib.request.Request(url, headers={"User-Agent": "orphan-domains/1D", "Accept": "text/x-fasta"})
    for _ in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
                lines = [ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">")]
                return "".join(lines)
        except urllib.error.URLError:
            time.sleep(1.0)
        except Exception:
            break
    return None

def compare_sequences(seq_struct: Optional[int], seq_up: Optional[str]) -> Tuple[Optional[bool], Optional[int], Optional[int], Optional[float]]:
    if seq_struct is None:
        return None, None, len(seq_up) if seq_up else None, None
    if not seq_up:
        return None, seq_struct, None, None
    len_up = len(seq_up)
    match = (seq_struct == len_up)
    ratio = (min(seq_struct, len_up) / max(seq_struct, len_up)) if max(seq_struct, len_up) > 0 else 1.0
    return match, seq_struct, len_up, ratio

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/afdb/index/swissprot_manifest.tsv")
    ap.add_argument("--acc-list", default="data/afdb/index/swissprot_accessions.txt")
    ap.add_argument("--pae-dir", default="data/afdb/pae")
    ap.add_argument("--out-dir", default="data/outputs/integrity")
    ap.add_argument("--n-sample", type=int, default=200)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    manifest = Path(args.manifest)
    acc_list = Path(args.acc_list)
    pae_dir = Path(args.pae_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    accs = read_accessions(acc_list)
    random.seed(args.seed)
    sample = accs if len(accs) <= args.n_sample else random.sample(accs, args.n_sample)

    summary_path = out_dir / "sample_verification.tsv"
    anomalies_path = out_dir / "anomalies.tsv"
    with summary_path.open("w") as s_out, anomalies_path.open("w") as a_out:
        header = [
            "uniprot","cif_gz","parse_ok","residues",
            "bfactor_min","bfactor_med","bfactor_max","plddt_in_0_100",
            "pae_present","pae_dim","pae_square_ok","uniprot_len",
            "seq_len_match","len_agreement_ratio","notes"
        ]
        s_out.write("\t".join(header) + "\n")
        a_out.write("\t".join(header) + "\n")

        for acc in sample:
            cif = find_cif_gz_for(acc, manifest)
            if cif is None:
                row = [acc,"","0","","","","","", "0","", "0","", "", "", "missing_manifest_entry"]
                s_out.write("\t".join(map(str,row)) + "\n")
                a_out.write("\t".join(map(str,row)) + "\n")
                log.warning("missing_manifest_entry", extra={"acc": acc})
                continue

            structure, resid_count, b_factors, msg = parse_structure_from_gz(cif)
            parse_ok = int(structure is not None and resid_count is not None and msg == "ok")
            bmin = f"{min(b_factors):.2f}" if b_factors else ""
            bmed = f"{statistics.median(b_factors):.2f}" if b_factors else ""
            bmax = f"{max(b_factors):.2f}" if b_factors else ""
            plddt_ok = int(bool(b_factors) and min(b_factors) >= 0.0 and max(b_factors) <= 100.0)

            pae_json = pae_dir / f"{acc}.json"
            pae_present, pae_dim, pae_msg = load_pae_json(pae_json)
            pae_square_ok = int(bool(pae_dim))

            up_seq = fetch_uniprot_fasta(acc)
            seq_match, len_struct, len_up, ratio = compare_sequences(resid_count, up_seq)
            seq_match_str = "" if seq_match is None else int(seq_match)
            ratio_str = "" if ratio is None else f"{ratio:.3f}"

            note_parts = []
            if msg != "ok": note_parts.append(msg)
            if pae_present and pae_msg != "ok": note_parts.append(pae_msg)
            if not pae_present: note_parts.append("pae_missing")
            # Flag size mismatch between PAE and structure (tolerate small gaps)
            if pae_dim and resid_count and abs(pae_dim - resid_count) > max(5, int(0.02 * max(pae_dim, resid_count))):
                note_parts.append("pae_struct_len_mismatch")
            if seq_match is False: note_parts.append("length_mismatch")
            notes = ",".join(note_parts) if note_parts else "ok"

            row = [
                acc, str(cif), str(parse_ok), str(len_struct or ""),
                bmin, bmed, bmax, str(plddt_ok),
                str(int(pae_present)), str(pae_dim or ""), str(pae_square_ok),
                str(len_up or ""), str(seq_match_str), str(ratio_str), notes
            ]
            s_out.write("\t".join(row) + "\n")
            if notes != "ok":
                a_out.write("\t".join(row) + "\n")

    log.info("verification_complete", extra={"summary": str(summary_path), "anomalies": str(anomalies_path)})
    print(f"Wrote {summary_path}")
    print(f"Wrote {anomalies_path}")

if __name__ == "__main__":
    main()
