# Biobakery Workflows in Snakemake

## Goals

This pipeline was built to have more fine control over the optimization of the components of Biobakery Workflows (how many cores/memory different parts used, when intermediate files were removed for disk space, etc). Snakemake makes these types of personalizations within different environments much easier compared to the ANADAMA manager.

In testing, a 60GB dataset utilizes around 500GB of disk space at most during the life of the job. Running the same dataset in the ANADAMA workflow with the same parameters resulted in ~1.5TB of disk space utilized.

## Installation

The first step is to download this repository to your local environment.

```
git clone https://github.com/ktmbiome-niaid/biobakery_snakemake.git
```

An `environment.yml` file was included in the repository for easy installation into the local environment.

With conda, create the new environment and install packages:

```
conda env create -n biobakery_snakemake -f environment.yml
```

The environment includes all dependencies necessary to run the pipeline, including snakemake.

Once the environment has been created, [there is one bug in kneaddata that we need to pre-empt](https://forum.biobakery.org/t/kneaddata-installed-with-conda-is-not-available/4147/1). To do this:

```
cd $CONDA_ENV_LOCATION/kneaddata/bin
ln -s ../share/trimmomatic/trimmomatic.jar .
```

There is a possibility that this file already exists in the intended directory, but we have encountered situations where it is not -- running these two lines of code confirms everything is in the right place.

### Databases

Next, we install the databases that biobakery uses throughout the pipeline. This process can take awhile, so we recommend submitting a job to the scheduler of your HPC. The commands for a SLURM batch script to download the databases used by the current workflow are below:

```
#!/usr/bin/bash
#SBATCH --mem=30g
#SBATCH --time=24:00:00

## KNEADDATA

conda activate biobakery_snakemake

kneaddata_database --download human_genome bowtie2 kneaddata_db
## This command appropriately formats the kneaddata database.

## METAPHLAN

metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db metaphlan_db

## HUMANN

humann_databases --download chocophlan full humann_db
humann_databases --download uniref uniref90_diamond humann_db
humann_databases --download utility_mapping full humann_db
```

## Setup Prior to Running the Pipeline

### Sequence Data

FASTQ files of single-end or paired-end reads should be combined into the same directory. Note the location of the directory so that you can add it to your config file in a few steps.

### Mapping File

To connect files in your directory of FASTQ files, this pipeline utilizes a "mapping file". This file should be formatted as follows so that it can be read into the workflow:

|#SampleID|ForwardFastqFile|ReverseFastqFile|
|---|---|---|
|GItract|S1_GItract_R1.fastq.gz|S1_GItract_R2.fastq.gz|
|Oral|S2_Oral_R1.fastq.gz|S2_Oral_R2.fastq.gz|
|Air|S3_Air_R1.fastq.gz|S3_Air_R2.fastq.gz|

The value in the `#SampleID` column will be used for subsequent naming of samples and files in processed results. We therefore recommend that this column not include special characters other than underscores.

### Config File

The config file provides instructions and unique locations to the Snakemake pipeline so that it can run successfully. For instance, this is where you will tell the pipeline where you saved the sequence data.

```
scratch_dir: "/lscratch/"
output_dir: "outputs/"

install_dir: "/my/installation/directory/"
map_file: "/my/installation/directory/mapping_file.txt"
input_dir: "/my/reads/directory/"
databaseDIR: "/directory/of/databases/"

paired_end: "False" # Needs to be string, looking for "True"
uniref_lev: "uniref90"
uniref_grp: "level4ec" # can choose rxn, go, ko, level4ec, pfam, eggnog
name: "ec"
units: "relab"

strainphlan: False

kd_db_loc: "kneaddata_db/"
kd_index: "hg37dec_v0.1"
mph_db_loc: "metaphlan_db/"
mph_index: "mpa_vJun23_CHOCOPhlAnSGB_202307"
```

Below provides a description of each of the variables:

- `scratch_dir`: Indicates the location of a directory where read/write will occur, typically the scratch directory on most servers.
- `output_dir`: The name of the directory where all results from the pipeline will be saved.
- `install_dir`: Directory where the biobakery_snakemake repository is saved.
- `map_file`: Location and name of the mapping file being used for analysis.
- `input_dir`: As noted earlier, location of the FASTQ reads for analysis.
- `databaseDIR`: Location of databases
- `paired_end`: Whether the pipeline is intended to run as single-end or paired-end. Choose `False` for single-end and `True` for paired-end.
- `uniref_lev`: Choose "uniref90" or "uniref50"
- `uniref_grp`: How HUMAnN groups uniref values; Choose "rxn", "go", "ko", "level4ec", "pfam", "eggnog"
- `name`: Linker to provide names for the new groups; Options are: "go", "ko", "ec", "pfam", "eggnog", matching `uniref_grp` order.
- `units`: Choose "relab" or "cpm"
- `strainphlan`: Whether to attempt StrainPhlAn
- `kd_db_loc`: Name chosen for the KneadData database
- `kd_index`: Name of the KneadData human genome index
- `mph_db_loc`: Name chosen for the MetaPhlAn Database
- `mph_index`: Version of the installed MetaPhlAn Database

## Modify the SBATCH script

Now that databases are installed, samples have been located, a mapping file has been created, and the config file has been set up for your analysis, let's edit the SLURM batch script. An example sbatch script has been provided (`example_sbatch_script.sh`). There are two key changes needed:

1. **SLURM resource requests.**

The header contains SLURM resource requests for the entire pipeline. Change or remove these requests as needed/required.

--time: How long, in hh:mm:ss, the pipeline should run.
--mem: The amount of memory that will be allocated to the entire job. This can be entirely dependent on the data being analyzed. 192GB is sufficient for most datasets. Some datasets may be able to go as low as 100GB.
--cpus-per-task: The number of CPUs to allocate. Computationally-intensive rules require a maximum of between 4-8 CPUs per task and sample, chosen to optimize the number of CPUs running at any given point in the pipeline.
--gres=lscratch: Optional. The HPC environment used for testing the pipeline requires scratch space requests.

2. **Config file location.**

In the script, you will see a line containing `export CONFIG_FILE="config.yaml"`. Provide the location or name of your config file.

You likely will not need to change this, but if you have moved the config.yaml file OR you have changed the name of the config.yaml (for instance, in the case where you want to run your samples as single-end in one job and paired-end in another job, creating two configs: `config_se.yaml` and `config_pe.yaml`).

### Submit the SBATCH script!
