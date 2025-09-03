
afdb_swissprot_cif: ## Download AFDB Swiss-Prot mmCIF tar (resume-safe)
	@URL=https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/swissprot_cif_v4.tar; \
	OUT=$(AFDB_DIR)/tar/swissprot_cif_v4.tar; \
	if command -v aria2c >/dev/null 2>&1; then aria2c -x8 -s8 -k1M -c -o $$(basename $$OUT) -d $$(dirname $$OUT) $$URL; \
	else wget -c -O $$OUT $$URL; fi

afdb_extract_cif: ## Extract .cif.gz files into data/afdb/cif_gz
	tar -xf $(AFDB_DIR)/tar/swissprot_cif_v4.tar -C $(AFDB_DIR)/cif_gz

afdb_index: ## Build accession manifest from filenames
	$(PY) - <<'PY'
from pathlib import Path
import re
cif_dir = Path("data/afdb/cif_gz")
idx_dir = Path("data/afdb/index"); idx_dir.mkdir(parents=True, exist_ok=True)
manifest = idx_dir / "swissprot_manifest.tsv"; accs = idx_dir / "swissprot_accessions.txt"
pat = re.compile(r"^AF-(?P<acc>[^-]+)-F\\d+-model_v4\\.cif\\.gz$")
rows, acc_list = [], []
for p in sorted(cif_dir.glob("*.cif.gz")):
    m = pat.match(p.name)
    if m:
        acc = m.group("acc"); rows.append((acc, str(p))); acc_list.append(acc)
manifest.write_text("uniprot\tcif_gz_path\n" + "\n".join(f"{a}\t{p}" for a,p in rows) + "\n")
accs.write_text("\n".join(acc_list) + "\n")
print(f"Wrote {manifest} and {accs} ({len(acc_list)}).")
PY

afdb_pae: ## Download PAE JSON for all accessions (parallel)
	./scripts/fetch_afdb_pae.sh

pfam_download: ## Fetch Pfam-A HMM and index
	cd data/pfam && wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && gunzip -f Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
