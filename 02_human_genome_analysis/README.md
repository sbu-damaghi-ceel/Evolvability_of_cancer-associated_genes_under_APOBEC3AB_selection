# CDUR analysis of human genome

## Requirement
For this analysis, you need to install CDUR previously developed in [Shapiro, Meier & MacCarthy 2018](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2161-y). You can found in [CDUR GitLab](https://gitlab.com/maccarthyslab/CDUR) page. If you followed the instruction but failed to run, follow [**How to install CDUR**](##-How-to-install-CDUR) instructions below

| Package | Version used in this study|
| --- | --- |
| Python | 3.6.13 (for CDUR), 3.9.16 |
| Biopython | 1.77 (for CDUR), 1.79 |
| CDUR | |

## Run CDUR
This will run CDUR for all sequences in `gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa`. Three replicates will be generated in separate folders.
```console
./run_cdur_refseq.sh
```

## Extract motif under-representation and mutational susceptibility values
This will generate `cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.csv` file containing motif under-representation and mutational susceptibility to draw figure.
```console
python3 extract_values_refseq_TC.py
```

## Draw Figure1a
One can use R script `draw_cdur_plot_refseq.r` to generate CDUR plot **Figure 1a** in the manuscript.

## How to install CDUR
CDUR was developed before `Biopython 1.78`, and utilize `Bio.Alphabet` module now removed.
Here we provide a way to install CDUR in `conda` environment.

### Make `conda` envirionment
```console
conda create --name cdur python=3.6
conda activate cdur
```

### Clone CDUR
```console
git clone https://gitlab.com/maccarthyslab/CDUR.git
```

### Install required package
```console
conda install -c conda-forge gsl --yes
conda install -c conda-forge numpy --yes
conda install -c conda-forge scipy --yes
conda install -c conda-forge pandas --yes
conda install -c conda-forge biopython=1.77 –yes
```

### Change makefile
Add `-I$${HOME}/miniconda3/envs/cdur/include -L$${HOME}/miniconda3/envs/cdur/lib` before `-o` option.

**Original**
```Makefile
all:
	g++ -O3 HotspotStatisticsReporter.cpp util.cpp StateVector.cpp Sequence.cpp Mutation.cpp scharff_utils.cpp ContinuousHistogram.cpp GeneticCode.cpp DiscreteHistogram.cpp MotifReference.cpp MotifIdentifier.cpp Motif.cpp EgnetProperties.cpp FrequencyDependentRandomizer.cpp Properties.cpp SequenceDataset.cpp StatsSampler.cpp BivariateNormalConditional.cpp generate.cpp analysis.cpp MotifMutationPair.cpp MotifMutationFrequency.cpp StatsSampler2Vars.cpp SimpleFastaReader.cpp RandomizedIota.cpp -I./ -I./tnt -o shmsim -lgsl -lgslcblas -lm -O

```
**Modified**
```Makefile
all:
	g++ -O3 HotspotStatisticsReporter.cpp util.cpp StateVector.cpp Sequence.cpp Mutation.cpp scharff_utils.cpp ContinuousHistogram.cpp GeneticCode.cpp DiscreteHistogram.cpp MotifReference.cpp MotifIdentifier.cpp Motif.cpp EgnetProperties.cpp FrequencyDependentRandomizer.cpp Properties.cpp SequenceDataset.cpp StatsSampler.cpp BivariateNormalConditional.cpp generate.cpp analysis.cpp MotifMutationPair.cpp MotifMutationFrequency.cpp StatsSampler2Vars.cpp SimpleFastaReader.cpp RandomizedIota.cpp -I./ -I./tnt -I$${HOME}/miniconda3/envs/cdur/include -L$${HOME}/miniconda3/envs/cdur/lib -o shmsim -lgsl -lgslcblas -lm -O

```
### Modify CDUR.py
From line 362, please modify original code as following to use `-o` and `-d` option.

**Original**
```Python
    if args.motifs is not None:
        os.system("~/CDUR/shmsim "+args.out_folder + seq_name[1:-1]+'_'+args.random_type+'.fasta ' + args.motifs + ' > ' + args.out_folder + seq_name[1:-1] + '_' + args.random_type + 'results.txt')

        if args.delete:
            try:
                os.remove(args.out_folder + seq_name[1:-1]+'.fas')
            except Exception:
                pass

            try:
                os.remove(f"{args.out_folder}{seq_name[1:-1]}_{args.random_type}.fasta")
            except Exception:
                pass
    else:
	    os.system("~/CDUR/shmsim "+seq_name[1:-1]+'_'+args.random_type+'.fasta > ' + seq_name[1:-1] + '_' + args.random_type + 'results.txt')

```

**Modified**
```Python
    if args.motifs is not None:
        os.system("~/CDUR/shmsim "+args.out_folder + seq_name[1:-1]+'_'+args.random_type+'.fasta ' + args.motifs + ' > ' + args.out_folder + seq_name[1:-1] + '_' + args.random_type + 'results.txt')

        if args.delete:
            try:
                os.remove(args.out_folder + seq_name[1:-1]+'.fas')
            except Exception:
                pass

            try:
                os.remove(f"{args.out_folder}{seq_name[1:-1]}_{args.random_type}.fasta")
            except Exception:
                pass
    else:
        os.system("~/CDUR/shmsim "+args.out_folder + seq_name[1:-1]+'_'+args.random_type+'.fasta > ' + args.out_folder + seq_name[1:-1] + '_' + args.random_type + 'results.txt')

        if args.delete:
            try:
                os.remove(args.out_folder + seq_name[1:-1]+'.fas')
            except Exception:
                pass

            try:
                os.remove(f"{args.out_folder}{seq_name[1:-1]}_{args.random_type}.fasta")
            except Exception:
                pass
```

Run following command in terminal to set `shmsim` location.

```console
PATH_TO_CDUR=/path/to/cdur/you/have/
sed -i "s+~/CDUR/shmsim+${PATH_TO_CDUR}/shmsim+g" ${PATH_TO_CDUR}/CDUR.py
```

## Run CDUR to all transcripts
```console
./run_cdur_all_gencode.sh
```
