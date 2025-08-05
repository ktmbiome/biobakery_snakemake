## Opening Statement
In the current implementation of our BioBakery pipeline, we utilized BioBakery Workflows, developed by the Huttenhower lab. The pipeline processes are managed by a tool called ANADAMA, which is also developed by the Huttenhower lab, though less-actively now. This means that we as the Nephele team have very limited ability to modify how components of the pipeline run, how files are manageed throughout, and how it interacts with the EC2. This has required that we use an excessively expensive EC2 instance type to be able to manage a normal amount of metagenomic data.

Therefore, the pipeline was rebuilt with all of the components using Snakemake as the pipeline manager. This means that we now have more fine control over the pipeline processes (ie, StrainPhlAn didn't *actually* work on most datasets in the BioBakery Workflows),what a user can provide to the pipeline, and we also have control over how files should be managed on the backend. We can remove unhelpful intermediate files to help save space, allowing us to spend computational disk space on the files that are required.

In addition, the report provided by the pipeline at the end has been re-designed into a quarto document, and now provides interactive outputs and plots.

Because of the total rewrite, the Snakefile is not properly integrated into the Nephele ecosystem and will require engineer expertise to integrate, including development of a main.py file. I am, however, available to answer any questions you may have!
## Software Tools, Scripts, Databases

### Tools
All tools that should be needed are provided in an [environment.yml](https://github.com/ktmbiome/biobakery_snakemake/blob/main/environment.yml) file. This pulls all the needed packages, including everything for the new and improved quarto document.

Channels need to include `biobakery` so that BioBakery-specific tools are downloaded.

Dependencies:
- python 3.9
- kneaddata 0.12.0
- metaphlan 4.1.0
- humann 3.9
- snakemake
- quarto (for interactive documents)
- r-kableextra
- r-heatmaply
- r-plotly
### Scripts
In this section, I will go through the scripts that are present in the github.

- [README.md](https://github.com/ktmbiome/biobakery_snakemake/blob/main/README.md) - this file is the readme for the repo, and provides technical information about downloading databases and setting up a run locally
- [Snakefile](https://github.com/ktmbiome/biobakery_snakemake/blob/main/Snakefile) - this is the main set of scripts. I'll go through each rule individually below.
- [config.yaml](https://github.com/ktmbiome/biobakery_snakemake/blob/main/config.yaml) - I created a temporary config file for someone to use.
- [environment.yml](https://github.com/ktmbiome/biobakery_snakemake/blob/main/environment.yml) - A file to define all of the software needed to run the pipeline
- [example_sbatch_script.sh](https://github.com/ktmbiome/biobakery_snakemake/blob/main/example_sbatch_script.sh) - Not useful for implementation, but I wrote a bash script to submit the snakemake job to a cluster.
- [utils.smk](https://github.com/ktmbiome/biobakery_snakemake/blob/main/utils.smk) - A utilities script to define a few necessary functions needed in the pipeline. This includes the `make_biom_from_tsv` script needed to convert metaphlan output to a BIOM file that is readable and acceptable by MicrobiomeDB and BIOM itself.
- [biobakery_scripts/count_features.py](https://github.com/ktmbiome/biobakery_snakemake/blob/main/biobakery_scripts/count_features.py) - A short python script from BioBakery Workflows that does feature counting. 
- [biobakery_scripts/get_counts_from_humann_logs.py](https://github.com/ktmbiome/biobakery_snakemake/blob/main/biobakery_scripts/get_counts_from_humann_logs.py) - Another short python script from BioBakery Workflows that parses the log to obtain the values used in the visualization pipeline tables.
- [biobakery_scripts/visualizations.qmd](https://github.com/ktmbiome/biobakery_snakemake/blob/main/biobakery_scripts/visualizations.qmd) - Quarto script written in R by Katie, heavily influenced by the BioBakery Workflows, to generate the visualizations HTML file.
### Databases
Documentation for downloading necessary databases is provided in the [README](https://github.com/ktmbiome/biobakery_snakemake?tab=readme-ov-file#databases) for the original Snakemake. In brief, after creating an environment that contains the necessary software, a `kneaddata`, `metaphlan` and `humann` database will need to be downloaded.
## Computational Requirements
For most datasets, this implementation of BioBakery should require an instance type more consistent with WGSA2. In testing, a 60GB dataset utilizes around 500GB of disk space at most during the life of the job. Running the same dataset in the ANADAMA workflow with the same parameters resulted in ~1.5TB of disk space utilized.
## Front-End Changes
This pipeline allows the user to make more selections than were previously possible. Below are the front-end selections that should be available to the user.
### BioBakery Card
Maintain the `Single-End` and `Paired-End` options.

Summary text should read: The bioBakery WGS pipeline utilizes the BioBakery suite of short read metagenomics tools including kneaddata, MetaPhlAn and HUMAnN to obtain assembly-free taxonomic and functional profiles.
### Pipeline Parameters:

| Option title     | Required? | UI type    | options (default in bold)                                                                                                                                          | Help-text if needed                                                                                                                                                                        |
| ---------------- | --------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Job Description  | Yes       | Short Text | -                                                                                                                                                                  | {none}                                                                                                                                                                                     |
| UniRef Level     | Yes       | Drop-down  | **"UniRef90"("uniref90")** OR "UniRef50"("uniref50")                                                                                                               | UniRef clusters proteins to reduce redundancy, either at 90% identity or 50% identity, respectively. We recommend UniRef90, but choose UniRef50 if you want broader functional annotation. |
| UniRef Group (+) | Yes       | Drop-down  | **"Enzyme Commission (EC)" ("level4ec")**, OR "Reaction" ("rxn") OR "Gene Ontology (GO)" ("go") OR "KEGG Orthology" ("ko") OR "Pfam" ("pfam") OR eggNOG ("eggnog") | Choose how you want HUMAnN3 to group reactions using commonly used databases.                                                                                                              |
| Units            | Yes       | Drop-down  | **"Relative Abundance" ("relab")** OR "Counts per Million" ("cpm")                                                                                                 | HUMAnN3 automatically normalizes the table before sharing the results. Indicate whether you would like to receive relative abundances or copies per million.                               |
| Run StrainPhlAn  | No        | Checkbox   | -                                                                                                                                                                  | Run StrainPhlAn to identify strains in metagenomic data                                                                                                                                    |
(+) This is a slightly difficult parameter. We will also need to link the "name" parameter in my config file with this, which aren't always the same. Some are easy, like uniref_grp: "go" matches to name: "go", but "rxn" matches with "", and "level4ec" matches to just "ec". But we need to make sure that this selection also fills in the correct name parameter.
## Snakemake Workflow
