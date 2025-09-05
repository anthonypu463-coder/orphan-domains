
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

orphan_qc: ## QC orphan candidates (coverage < 0.20 + pLDDT stats; structural subset)
	$(PY) scripts/qc_orphan_candidates.py

orphan_datasets: ## Produce canonical <20% set + strict (~26k) + strict-structural subsets
	$(PY) scripts/compile_orphan_datasets.py

orphan_zero_struct70: ## Build 0% Pfam & mean pLDDT ≥70 primary dataset
	$(PY) scripts/make_orphan_zero_struct70.py

cluster: ## MMseqs2 clustering of orphan_zero_struct70 (40% ID, cov 0.8 of shorter)
	scripts/cluster_orphanZ70.sh

cluster26k: ## Cluster the ~26k set at 40% ID (single-linkage, cov 0.8 of shorter)
	scripts/cluster_strict26k.sh

cluster26k_reps: ## Select cluster representatives by highest mean pLDDT (tie: length)
	$(PY) scripts/select_cluster_representatives.py

cluster26k_check: ## Analyze 26k cluster sizes; flag large clusters (>100); write histogram & summary
	$(PY) scripts/check_cluster_sizes.py

cluster26k_outputs: ## Build member→cluster assignments + size summary (TSV/JSON)
	$(PY) scripts/compile_cluster_outputs.py

stage_structures: ## Stage AlphaFold structures (symlinks) for analysis; SET=reps_113k|reps_26k|all_113k, LIMIT=N optional
	$(PY) scripts/prepare_structure_workspace.py --set $${SET:-reps_113k} ${LIMIT:+--limit $${LIMIT}}

struct_consistency: ## Check sequence↔structure consistency; SET=reps_113k|reps_26k|all_113k, LIMIT=N optional
	$(PY) scripts/check_sequence_structure_consistency.py --set $${SET:-reps_113k} ${LIMIT:+--limit $${LIMIT}}

struct_confidence: ## Verify pLDDT (B-factor) and ensure PAE JSON exists; SET=orphanZ70_reps|strict26k_reps|orphanZ70_all
	$(PY) scripts/collect_confidence_data.py --set $${SET:-orphanZ70_reps}

struct_features: ## Extract per-residue features (Cα coords + pLDDT) to NPZ; SET=orphanZ70_reps|strict26k_reps|orphanZ70_all
	$(PY) scripts/extract_structure_features.py --set $${SET:-orphanZ70_reps} ${LIMIT:+--limit $${LIMIT}} ${OVERWRITE:+--overwrite}

segment: ## 7A: find low-confidence pLDDT regions (pLDDT < THR, default 50); SET=orphanZ70_reps|strict26k_reps|orphanZ70_all
	$(PY) scripts/find_low_plddt_regions.py --set $${SET:-orphanZ70_reps} --threshold $${THR:-50} --min-run $${MINRUN:-3} ${LIMIT:+--limit $${LIMIT}}

plddt_scan: ## 7A — find low-confidence spans (pLDDT < 50); SET=orphanZ70_reps|strict26k_reps|orphanZ70_all
	$(PY) scripts/find_low_plddt_regions.py --set $${SET:-orphanZ70_reps} ${LIMIT:+--limit $${LIMIT}}

segment_plddt: ## 7B — segment by pLDDT cuts (midpoints), keep fragments >=50 aa; SET=...
	$(PY) scripts/segment_by_plddt_cuts.py --set $${SET:-orphanZ70_reps} ${LIMIT:+--limit $${LIMIT}}

segment_refine: ## 7C — enforce min fragment length (>=50 aa); merge interior shorts, drop terminal tails
	$(PY) scripts/enforce_min_fragment_length.py --set $${SET:-orphanZ70_reps} ${LIMIT:+--limit $${LIMIT}} ${MIN:+--min-fragment $${MIN}}
