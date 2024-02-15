# Compare cancer and control genes

## Download COSMIC cancer census
We download a list of cancer-associated genes from [COSMIC Data Downloads](https://cancer.sanger.ac.uk/cosmic/archive-download) page. One can download directly from the page or use script to download the file.
```console
./download_cosmic_cancer_gene_census.sh
```

## Annotate gene class to CDUR result file
We annotate gene class (cancer, control, or regular) to CDUR result file and generate `cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.geneclass.csv` to draw comparison figure.
```console
python3 annotate_cancer_genes.py
```

## Draw comparison figure
One can use R script `cancer_vs_control.r` to draw **Figure 2a**, and also perform statistical tests.