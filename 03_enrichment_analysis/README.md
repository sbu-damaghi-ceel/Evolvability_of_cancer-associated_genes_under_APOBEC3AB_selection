# Enrichment analysis

## Extract the genes from CDUR plot
One can use R script `gene_enrichment_analysis_refseq.r` to extract top-left and bottom-right genes of **Figure 1a**. This script will generate `top_left.txt` and `bottom_right.txt`.

## ShinyGO
We used [ShinyGO 0.77](http://bioinformatics.sdstate.edu/go/) for gene enrichment analysis. We used default parameters as in below table. We downloaded the result `Top Pathways shown above`.
| Parameter | Value |
| --- | --- |
| Species | Human |
| FDR cutoff | 0.05 |
| \# pathways to show | 20 |
| Pathway size: Min. | 2 |
| Pathway size: Max. | 2000 |
| Remove redundancy | checked |
| Abbreviate pathways | checked |

## Draw Figure 1b and 1c
One can use R script `gene_enrichment_figures.R` to draw **Figure 1b** and **1c**. Gene enrichment analysis results in supplementary figures were directly downloaded from ShinyGO result page.