# AE1712-viromes
Analysis of the AE1712 Long Read virome dataset


# Analysis for Figures

## Supplementary Figure XX (new) 
* Input file: `data/alignment_shortread_contigs_gr1kb_to_all_LRR-sorted-by-target-name.txt`
* Parsing file: `scripts/calc_breakages.py`
* Run parameters:
```
python scripts/calc_breakages.py    --alignment_file data/alignment_shortread_contigs_gr1kb_to_all_LRR-sorted-by-target-name.txt \
                                    --out_file output/breakages.csv
                                    --min_pct_align 70
```
The aim of this analysis is to understand whether or not viruses that are assembled through long reads are in the short read fraction, just in fragmented components.

All short read contigs >1kb were aligned against viral population representatives from populations where the longest member (the representative) was derived from a long-read assembly. The `PAF` output from `minimap2` was parsed with the script `calc_breakages.py` to identify regions of the representative that were not covered by short read assemblies. 


