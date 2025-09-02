# Orphan Domains (Orphan Proteins)

Reproducible pipeline to discover novel protein domains from AFDB Swiss-Prot.

## Baseline environment (pinned)
- Python 3.11 (conda-forge)
- MMseqs2 v14
- HMMER 3.x
- Foldseek v1 (2023)
- DPAM v2.0

## Setup
mamba env create -f environment.yml
conda activate orphan-domains
python -m orphan.cli env_check --log-format text | tee docs/baseline_env.txt

## Layout
data/ (gitignored), docs/, src/orphan/, scripts/, tests/
