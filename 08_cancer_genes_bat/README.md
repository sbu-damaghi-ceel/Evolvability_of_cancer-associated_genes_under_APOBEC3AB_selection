# Bat genome cancer vs control comparison

## Download *Pteropus alecto* equivalent cancer-associated and control genes
```console
./make_list_orthologs.py
``` 

## Annotate gene class to bat CDUR result file
```console
python3 bat_annotate_cancer_genes.py
```

## Draw Figure 3b
One can use R script `bat_cancer_vs_control.r` to generate CDUR plot **Figure 3b** in the manuscript.