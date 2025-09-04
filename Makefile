
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

pfam_sequences: ## Build FASTA from AFDB mmCIFs
	$(PY) scripts/extract_sequences.py --workers 8

pfam_split: ## Split FASTA into shards
	scripts/split_fasta.sh

pfam_scan: ## Run hmmscan with Pfam-A --cut_ga across shards (parallel)
	JOBS?=8 THREADS_PER?=1 scripts/run_hmmscan.sh

pfam_parse: ## Parse domtblout into TSV and summarize coverage
	$(PY) scripts/parse_hmmscan_domtbl.py
	$(PY) scripts/summarize_pfam_coverage.py

pfam_resolve: ## Resolve overlaps: non-overlapping Pfam domains per protein
	$(PY) scripts/resolve_domain_overlaps.py

pfam_coverage: ## Compute Pfam coverage from non-overlapping domains
	$(PY) scripts/compute_pfam_coverage.py

pfam_finalize: ## Final non-overlapping Pfam architecture (optional same-clan trimming)
	$(PY) scripts/finalize_domain_architecture.py --trim-policy same-clan --min-fragment 30

pfam_arch_strings: ## Compile domain architecture strings + unannotated segments
	$(PY) scripts/compile_architecture_strings.py

pfam_arch_enhanced: ## Enhanced final annotations (domains + gap coordinates)
	$(PY) scripts/compile_domain_annotations_final.py

orphan_select: ## Select orphan-rich proteins (coverage < 0.20)
	$(PY) scripts/select_orphan_candidates.py

orphan_compile: ## Build orphan sequences FASTA + metadata from candidates
	$(PY) scripts/compile_orphan_sets.py
