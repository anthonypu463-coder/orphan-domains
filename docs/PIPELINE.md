# Pipeline — Steps 1–8 frozen
- 1–3: Pfam scan → non-overlapping Pfam map
- 4: Orphan selection (<20% coverage); subset defaults fixed
- 5: Clustering (40% ID)
- 6: Structures staged; seq↔struct OK; pLDDT/PAE present; features extracted
- 7: pLDDT segmentation (cut <50), ≥50aa filter, edge overrides; PAE validate/merge; tail trimming
- 8: Final domains minted (OD IDs), logs and summaries written

Next (not in this commit): 9 DPAM (ambiguous only), 10 ≥95% de-dup, 11 novelty, 12–13 graphs & splits.
