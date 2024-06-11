import os
import pandas as pd
import csv
import biom
import numpy as np
from typing import Any, Callable, Dict, List

def read_mapping_file():
    samples = {}
    with open(config["map_file"]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['#SampleID'] == '':
                continue
            samples[row['#SampleID']] = {}
            if row.get('ForwardFastqFile'):
                samples[row['#SampleID']]["f"] = config["input_dir"] + row["ForwardFastqFile"]
            if row.get('ReverseFastqFile'):
                samples[row['#SampleID']]["r"] = config["input_dir"] + row["ReverseFastqFile"]
    return samples

filename = config["map_file"]
samples = read_mapping_file()
print(samples)

def get_fastq_gz(wildcards):
    return samples[wildcards.sample]

def get_var(bash_var):
    var = os.environ.get(bash_var)
    return int(var)

threads_per_job = get_var('SLURM_CPUS_PER_TASK')
print("threads available:", threads_per_job)

mem = os.environ.get('SLURM_MEM_PER_NODE')
if isinstance(mem, str) and len(mem) > 0:
    total_mem_gb = int(get_var('SLURM_MEM_PER_NODE') / 1024)
else:
    total_mem_gb = int(get_var('SLURM_MEM_PER_CPU') * threads_per_job  / 1024 )
print("memory available:", total_mem_gb, "Gb")

tmpdir = str(get_var('SLURM_JOB_ID'))
print("jobID:", tmpdir)

wDIR = os.getcwd()
print("working directory:", wDIR)


def make_biom_from_tsv(infile, outfile, mapfile=config["map_file"]):
    """
    Make a BIOM file from the tsv files created by biobakery (metaphlan and humann)

    Params
    -----
    outputs_dir = location of outputs (set by arguments)
    file_path = file dir of tsv file within outputs (set by biobakery)
    infile = input filename (present within outputs_dir + file_path)
    outfile = resulting BIOM file name (placed within outputs_dir + file_path)
    samples = samples for insertion into the BIOM file

    Returns
    -----
    BIOM-formatted file
    """

    def filter_taxonomy(line):
        if not ("clade_name" in line or "t__" in line):
            return
        return line.replace("|", ";")

    def filter_gene(line):
        if "|" in line:
            return
        return line.replace("_Abundance", "").replace("-RPKs", "")

    def filter_text_file(rawfile: List, filtering_function: Callable[[str], str]):
        for line in f:
            filtered_line = filtering_function(line)
            if not filtered_line:
                continue
            rawfile.append(filtered_line)

    rawfile = []
    with open(str(infile), "r", newline="\n") as f:
        first_word = f.readline().split("\t")[0].strip()
    is_taxonomy = first_word.startswith("#mpa")
    with open(str(infile), "r", newline="\n") as f:
        filter_text_file(rawfile, filter_taxonomy if is_taxonomy else filter_gene)
    data = list(zip(*(line.strip().split("\t") for line in rawfile)))
    data.append(data[0])
    data2 = np.array(data)
    count_data = np.transpose(data2[1:-1, 1:])
    observ_ids = list(data2[0][1:])
    sample_ids = list(np.transpose(data2)[0][1:-1])
    updt_first_word = "taxonomy" if is_taxonomy else first_word
    observ_dict = []
    for x in observ_ids:
        observ_dict.append({updt_first_word: x.split(";")})

    with open(str(mapfile), "r", newline="\n") as m_in:
        sample_dict = list(csv.DictReader(m_in, delimiter="\t"))
    sample_dict.sort(key=lambda a: a["#SampleID"])

    biom_tab = biom.Table(
        count_data,
        sample_ids=sample_ids,
        observation_ids=observ_ids,
        observation_metadata=observ_dict,
        sample_metadata=sample_dict
    )
    with open(str(outfile), "w") as f:
        biom_tab.to_json("Nephele", direct_io=f)
    print(f"Created BIOM file: {outfile}")
