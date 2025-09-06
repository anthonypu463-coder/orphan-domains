# Orphan Domains (AFDB Swiss-Prot)

A reproducible pipeline and dataset of **novel protein domains** from AlphaFold DB (Swiss-Prot).

**Current freeze:** Steps **1–8** complete (final domains minted with stable IDs).  
Next (not yet included here): **Step 9 (DPAM)**, **Step 10 (≥95% de-dup)**, **Step 11 (novelty)**, **Steps 12–13 (graphs & ML splits)**.

## Locked thresholds
- pLDDT cut: **<50** → candidate linker
- Minimum fragment: **≥50 aa**
- PAE interface heuristic: **<5 Å** suggests merge, **>15 Å** supports split (window=10)
- All decisions logged in `data/outputs/segments/orphanZ70_reps/*.json`

## Key outputs now
- Final domains (OD IDs):  
  `data/outputs/domains/orphanZ70_reps/domains_final.{tsv,fasta,jsonl}`  
  `data/outputs/domains/orphanZ70_reps/domains_summary.json`
- Segmentation/PAE refinement logs:  
  `data/outputs/segments/orphanZ70_reps/segments_refined.tsv`  
  `data/outputs/segments/orphanZ70_reps/segments_trimmed.tsv`  
  `data/outputs/segments/orphanZ70_reps/segments_refine_summary.json`  
  `data/outputs/segments/orphanZ70_reps/segments_trim_summary.json`  
  `data/outputs/segments/orphanZ70_reps/pae_interface_summary.json`
- Registry: `data/outputs/registry/domain_registry.tsv`

