#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT"

pass=true
ok(){ printf "PASS  %s\n" "$*"; }
fail(){ printf "FAIL  %s\n" "$*\n"; pass=false; }

# A) Repo & baseline files
[ -d .git ] && ok "git repo exists" || fail "git repo missing"
branch=$(git symbolic-ref --quiet --short HEAD || echo "?")
[ "$branch" = "main" ] && ok "on main branch" || fail "not on main (got: $branch)"
[ -s README.md ] && ok "README present" || fail "README missing"
[ -s CITATION.cff ] && ok "CITATION.cff present" || fail "CITATION.cff missing"
[ -s Makefile ] && ok "Makefile present" || fail "Makefile missing"
[ -s pyproject.toml ] && ok "pyproject.toml present" || fail "pyproject.toml missing"
[ -s src/orphan/cli.py ] && ok "CLI present" || fail "CLI missing"

# B) Dataset versioning + IDs
[ -s docs/dataset_versions.yaml ] && ok "dataset_versions.yaml present" || fail "dataset_versions.yaml missing"
[ -s data/outputs/registry/domain_registry.tsv ] && ok "domain_registry.tsv present" || fail "domain_registry.tsv missing"
[ -s data/outputs/registry/domain_id_counter.json ] && ok "domain_id_counter.json present" || fail "domain_id_counter.json missing"

# C) Data retrieval (AFDB + PAE + Pfam)
acc_file="data/afdb/index/swissprot_accessions.txt"
mani_file="data/afdb/index/swissprot_manifest.tsv"
[ -s "$acc_file" ] && ok "accession list present" || fail "accession list missing"
[ -s "$mani_file" ] && ok "manifest present" || fail "manifest missing"

acc_total=$(wc -l < "$acc_file" | tr -d '[:space:]')
cif_count=$(find data/afdb/cif_gz -maxdepth 1 -name '*.cif.gz' | wc -l | tr -d '[:space:]')
pae_count=$(find data/afdb/pae -maxdepth 1 -name '*.json' | wc -l | tr -d '[:space:]')

[[ "$cif_count" -gt 0 ]] && ok "CIF files extracted (~$cif_count)" || fail "no CIF files found"
[[ "$pae_count" -gt 0 ]] && ok "PAE JSON present (~$pae_count)" || fail "no PAE JSON found"

for f in data/pfam/Pfam-A.hmm data/pfam/Pfam-A.hmm.h3{m,i,f,p}; do
  [ -s "$f" ] || { fail "Pfam index missing: $(basename "$f")"; pfam_ok=false; }
done
${pfam_ok:=true} && ok "Pfam HMM pressed (h3m/h3i/h3f/h3p)"

# D) Integrity verification outputs
[ -s data/outputs/integrity/sample_verification.tsv ] && ok "integrity summary present" || fail "integrity summary missing"
anom="data/outputs/integrity/anomalies.tsv"
if [ -s "$anom" ]; then
  an_lines=$(( $(wc -l < "$anom") - 1 ))
  if (( an_lines == 0 )); then ok "anomalies: none in sample"; else fail "anomalies: $an_lines rows"; fi
else
  fail "anomalies.tsv missing"
fi

# E) Organization & reference manifest
[ -d data/afdb/by_acc ] && ok "by-acc directory exists" || fail "by-acc directory missing"
ref="data/afdb/index/reference_manifest.tsv"
if [ -s "$ref" ]; then
  rows=$(( $(wc -l < "$ref") - 1 ))
  if (( rows > 0 )); then ok "reference_manifest.tsv rows: $rows"; else fail "reference_manifest.tsv empty"; fi
else
  fail "reference_manifest.tsv missing"
fi

# Summary
echo "-----"
$pass && { echo "ALL CHECKS PASS"; exit 0; } || { echo "CHECKS FAILED"; exit 1; }
