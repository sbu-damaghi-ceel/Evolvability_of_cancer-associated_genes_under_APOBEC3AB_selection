# Data downloading and preprocessing

## Downloading GENCODE 
To download GENCODE protein coding transcripts sequence, change `/path/to/download/` in the file `download_genecode.sh` to the path you want to download tha file and run the script. 

```console
./download_genecode.sh
```

## Preprocessing GENCODE
The file download with above script contains not only coding sequence (CDS) but also 5' untranslated region (5'UTR) and 3' untranslated region (3'UTR). 
In addition, several transcripts have deprecated coding sequences. This scripts extract only CDS portion of valid transcripts that can produce polypeptide.

### Requirements
| Package | Version used in this study|
| --- | --- |
| Python | 3.9.16 |
| Biopython | 1.79 |

```console
python3 extract_cds_only.py
```

### Reducing number of transcripts
This step we used biomaRt package to get information from Ensembl, and use the information to reduce the number of transcripts to search.
We used two different criteria, 1) sequences marked as *canonical*, 2) sequences that have NCBI RefSeq ID. 

**You need to change `base_path` of all scripts file accordingly**

First, run `extract_name_id.sh` to generate CSV file containing only name and Ensembl transcript ID from `gencode.v40.pc_transcripts.nopary.cdsonly.fa`. 

```console
./extract_name_id.sh
```

Second, run `reduce_num_sequence.R` to retrieve Ensembl transcript ids. 

Lastly, run `reduce_num_transcripts.py` to generate FASTA file.

```console
python3 reduce_num_transcripts.py
```

This will generate two FASTA files, `gencode.v40.pc_transcripts.nopary.cdsonly.canonical.fa` and `gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa`