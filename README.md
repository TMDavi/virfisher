# virfisher

Virfisher is a snakemake pipeline for the exploration of viral sequences in metagenomics data. It performs the viral sequence identification, annotation, taxonomic assignment, deduplication, host assignment, life style prediction and relative abundance calculation

# Instalation

1. Clone the repository

```shell
git clone https://github.com/TMDavi/virfisher.git

```
2. Install Miniconda3

3. Install the virfisher environment

```shell
conda env create -f virfisher_env.yaml
conda activate virfisher_env
```

4. Install the tool databases

- Run the setup_databases.py script
```shell

virfisher set_databases --db_dir path/to/dir
```

# Easy run

1. You place your metagenomic reads into a single folder and run using the input folder parameter

```shell
python virfisher.py --input_folder data/ --outdir outtest/ --threads 5
```

2. Another way to run the pipeline is to use a csv file with the sample paths as an input

```shell
python virfisher.py --input_sheet example/samplesheet.csv --outdir outtest/ --threads 5
```

