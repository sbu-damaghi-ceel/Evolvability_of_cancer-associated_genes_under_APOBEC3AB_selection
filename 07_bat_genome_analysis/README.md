# Bat genome analysis

## Download a list of all genes of *Pteropus alecto*
We used [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (CLI v13.x API v1) to download a list of *Pteropus alecto* genes in NCBI. The shell script also run Python script `get_short_info_from_summary.py` to preprocess `.json` file and generate `.tsv` file. 
```console
./download_ncbi_palecto.sh
``` 

## Download valid transcripts' coding sequences
```console
python3 download_valid_seqs_getrefseq.py
```

## Run CDUR
```console
./run_cdur_refseq_palecto.sh
```

## Extract motif under-representation and mutational susceptibility values
```consoel
python3 extract_values_refseq_TC_palecto.py
```

## Draw Figure 3a
One can use R script `draw_cdur_plot_refseq_palecto.r` to generate CDUR plot **Figure 3a** in the manuscript.