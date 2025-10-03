configfile: "config.yaml"
include: "utils.smk"

TMPDIR = config["scratch_dir"] + tmpdir + "/"
WORKDIR = wDIR
KD_OUTPUTS=config["output_dir"] + "kneaddata/"
MPHL_OUTPUTS=config["output_dir"] + "metaphlan/"
HMN_OUTPUTS=config["output_dir"] + "humann/"
STRN_OUTPUTS=config["output_dir"] + "strainphlan/"
MERGED="merged/"
MAIN="main/"
COUNTS="counts/"
BIOMS="biom/"
RELEVS=config["units"] + "/"

pipeline_output = [
    TMPDIR + KD_OUTPUTS + MERGED + "kneaddata_read_count_table.tsv",
    TMPDIR + MPHL_OUTPUTS + MERGED + "metaphlan_taxonomic_profiles.tsv",
    expand(TMPDIR + HMN_OUTPUTS + MAIN + "{sample}_pathabundance.tsv", sample=samples),
    TMPDIR + HMN_OUTPUTS + MERGED + config["name"] + "s_" + config["units"] + ".tsv",
    TMPDIR + HMN_OUTPUTS + MERGED + config["name"] + "s_" + config["units"] + "_renamed.tsv",
    TMPDIR + MPHL_OUTPUTS + BIOMS + "metaphlan_taxonomic_profiles.biom",
    TMPDIR + HMN_OUTPUTS + BIOMS + "humann_pathabundance_" + config["units"] + ".biom",
    TMPDIR + HMN_OUTPUTS + BIOMS + "humann_" + config["name"] + "_" + config["units"] + ".biom",
    TMPDIR + HMN_OUTPUTS + MERGED + "humann_read_and_species_count_table.tsv",
    TMPDIR + HMN_OUTPUTS + COUNTS + "humann_feature_counts.tsv"
]

if len(samples) > 2:
    pipeline_output.extend([
        TMPDIR + config["output_dir"] + "visualizations.html"
    ])

if config["strainphlan"]:
    pipeline_output.extend([
        TMPDIR + STRN_OUTPUTS + "clade_names/print_clades_only.tsv"
    ])

rule all:
    input:
        pipeline_output

rule kneaddata:
    input:
        unpack(get_fastq_gz)
    output:
        reads=temp(TMPDIR + KD_OUTPUTS + MAIN + "{sample}.fastq"),
        log=TMPDIR + KD_OUTPUTS + MAIN + "logs/{sample}.log"
    params:
        database=config["databaseDIR"] + config["kd_db_loc"] + config["kd_index"],
        outdir=TMPDIR + KD_OUTPUTS + MAIN,
        logdir=TMPDIR + KD_OUTPUTS + MAIN + "logs/",
        log=TMPDIR + KD_OUTPUTS + MAIN + "logs/{sample}.log",
        paired=config["paired_end"]
    threads: 8
    resources:
        mem_mb=30000
    shell:
        """
        cmd="kneaddata "
        if [ {params.paired} == "True" ]; then
            cmd+="--input1 {input.f} --input2 {input.r} "
        else
            cmd+="--unpaired {input.f}"
        fi
        cmd+=" --log {params.log} \
        --reference-db {params.database} --max-memory {resources.mem_mb}m \
        --output {params.outdir} --run-fastqc-end --log-level INFO --cat-final-output \
        --remove-intermediate-output -t {threads} --output-prefix {wildcards.sample}"

        echo "$cmd"

        eval "$cmd"

        rm {params.outdir}/{wildcards.sample}_*.fastq
        """


rule metaphlan:
    input:
        reads=rules.kneaddata.output.reads
    output:
        bt2=temp(TMPDIR + MPHL_OUTPUTS + MAIN + "{sample}.bt2"),
        profile=TMPDIR + MPHL_OUTPUTS + MAIN + "{sample}_profile.tsv",
        sams = temp(TMPDIR + STRN_OUTPUTS + MAIN + "{sample}.sam.bz2")
    params:
        db=config["databaseDIR"] + config["mph_db_loc"],
        index=config["mph_index"]
    threads: 10
    shell:
        """
        metaphlan --nproc {threads} -s {output.sams} --bowtie2out {output.bt2} --input_type fastq \
        --bowtie2db {params.db} --index {params.index} -o {output.profile} \
        {input.reads}
        """

rule humann:
    input:
        reads=rules.kneaddata.output.reads,
        profile=rules.metaphlan.output.profile
    output:
        genefam = TMPDIR + HMN_OUTPUTS + MAIN + "{sample}_genefamilies.tsv",
        pathabun= TMPDIR + HMN_OUTPUTS + MAIN + "{sample}_pathabundance.tsv",
        pathcov = TMPDIR + HMN_OUTPUTS + MAIN + "{sample}_pathcoverage.tsv",
        log = TMPDIR + HMN_OUTPUTS + MAIN + "{sample}.log"
    params:
        outdir = TMPDIR + HMN_OUTPUTS + MAIN,
        db=config["databaseDIR"] + "humann_db/",
        index=config["mph_index"]
    threads: 16
    shell:
        """              
        humann -i {input.reads} -o {params.outdir} \
        --threads {threads} --o-log {output.log} --diamond-options "-b 2 -c 6" \
        --input-format fastq --search-mode uniref90 --minpath on \
        --nucleotide-database {params.db}chocophlan/ --protein-database {params.db}uniref/ \
        --taxonomic-profile {input.profile} --remove-temp-output --metaphlan-options "-x {params.index}"
        """

rule strainphlan_markers:
    input:
        rules.metaphlan.output.sams
    output:
        TMPDIR + STRN_OUTPUTS + "consensus_markers/{sample}.json.bz2"
    params:
        outdir= TMPDIR + STRN_OUTPUTS + "consensus_markers/",
        db=config["databaseDIR"] + config["mph_db_loc"] + "/" + config["mph_index"] + ".pkl",
    threads: 8 # can perhaps be scaled down? no evidence of getting scaled down
    shell:
        """
        mkdir -p {params.outdir}

        sample2markers.py -i {input} -o {params.outdir} --nproc {threads} \
            -d {params.db}
        """

rule strainphlan_clades:
    input:
        expand(TMPDIR + STRN_OUTPUTS + "consensus_markers/{sample}.json.bz2", sample=samples)
    output:
        TMPDIR + STRN_OUTPUTS + "clade_names/print_clades_only.tsv"
    params:
        db=config["databaseDIR"] + config["mph_db_loc"] + "/" + config["mph_index"] + ".pkl",
        dir=TMPDIR + STRN_OUTPUTS
    threads: 10
    shell:
        """
        mkdir -p {params.dir}/clade_names
        strainphlan -s {input} -o {params.dir}/clade_names -n {threads} --print_clades_only -d {params.db}
        if [[ $(wc -l < {params.dir}/clade_names/print_clades_only.tsv) -ge 2 ]]; then
            sed 1d {params.dir}/clade_names/print_clades_only.tsv | cut -f1 | while read CLADE; do 
                echo "Running clade $CLADE"
                mkdir -p {params.dir}/$CLADE
                strainphlan -s {params.dir}/consensus_markers/*.json.bz2 -o {params.dir}/$CLADE -n {threads} -c $CLADE -d {params.db}
            done
        else
            echo "No eligible clades were identifed (clades must be present at 80% coverage in at least 4 samples). StrainPhlAn will not run."
        fi
        """

# ------  Rules for combining files from previous rules ---------

rule kneaddata_read_count_table:
    input:
        files = expand(TMPDIR + KD_OUTPUTS + MAIN + "logs/{sample}.log", sample=samples)
    output:
        table=TMPDIR + KD_OUTPUTS + MERGED + "kneaddata_read_count_table.tsv"
    params:
        logdir = TMPDIR + KD_OUTPUTS + MAIN + "logs/"
    threads: 1
    shell:
        """
        kneaddata_read_count_table --input {params.logdir} --output {output.table}
        """

rule combine_metaphlan:
    input:
        expand(TMPDIR + MPHL_OUTPUTS + MAIN + "{sample}_profile.tsv", sample=samples)
    output:
        merged_file=TMPDIR + MPHL_OUTPUTS + MERGED + "metaphlan_taxonomic_profiles.tsv"
    threads: 1
    shell:
        """
        merge_metaphlan_tables.py {input} > {output.merged_file}
        """

rule humann_join_relev:
    input:
        expand(TMPDIR + HMN_OUTPUTS + MAIN + "{sample}_pathcoverage.tsv", sample=samples)
    params:
        dir = TMPDIR + HMN_OUTPUTS + MAIN,
        units= config["units"]
    output:
        comb_pathcov= TMPDIR + HMN_OUTPUTS + MERGED + "pathcoverage.tsv",
        comb_genefam= TMPDIR + HMN_OUTPUTS + MERGED + "genefamilies.tsv",
        comb_pathabd= TMPDIR + HMN_OUTPUTS + MERGED + "pathabundance.tsv",
        relev_pathcov= TMPDIR + HMN_OUTPUTS + RELEVS + "pathcoverage_" + config["units"] + ".tsv",
        relev_genefam= TMPDIR + HMN_OUTPUTS + RELEVS + "genefamilies_" + config["units"] + ".tsv",
        relev_pathabd= TMPDIR + HMN_OUTPUTS + RELEVS + "pathabundance_" + config["units"] + ".tsv"
    threads: 1
    shell:
        """
        humann_join_tables -i {params.dir} --file_name "pathcoverage" -o {output.comb_pathcov}
        humann_renorm_table -i {output.comb_pathcov} --units {params.units} -o {output.relev_pathcov}
        humann_join_tables -i {params.dir} --file_name "genefamilies" -o {output.comb_genefam}
        humann_renorm_table -i {output.comb_genefam} --units {params.units} -o {output.relev_genefam}
        humann_join_tables -i {params.dir} --file_name "pathabundance" -o {output.comb_pathabd}
        humann_renorm_table -i {output.comb_pathabd} --units {params.units} -o {output.relev_pathabd}
        """

rule humann_regroup:
    input:
        rules.humann_join_relev.output.relev_genefam,
    output:
        TMPDIR + HMN_OUTPUTS + MERGED + config["name"] + "s_" + config["units"] + ".tsv"
    params:
        group=expand("{levs}_{grps}", levs=config["uniref_lev"], grps=config["uniref_grp"]),
        file=config["databaseDIR"] + "humann_db/utility_mapping/map_" + config["uniref_grp"] + "_" + config["uniref_lev"] + ".txt.gz"
    threads: 1
    shell:
        """
        humann_regroup_table --input {input} -o {output} -c {params.file}
        """

rule humann_rename:
    input:
        TMPDIR + HMN_OUTPUTS + MERGED + config["name"] + "s_" + config["units"] + ".tsv"
    output:
        TMPDIR + HMN_OUTPUTS + MERGED + config["name"] + "s_" + config["units"] + "_renamed.tsv"
    params:
        name=config["name"],
        file=config["databaseDIR"] + "humann_db/utility_mapping/map_" + config["name"] + "_name.txt.gz"
    shell:
        """
        humann_rename_table --input {input} -o {output} -c {params.file}
        """

rule humann_counts_from_logs:
    input:
        expand(TMPDIR + HMN_OUTPUTS + MAIN + "{sample}.log", sample=samples)
    output:
        TMPDIR + HMN_OUTPUTS + MERGED + "humann_read_and_species_count_table.tsv"
    params:
        script_loc=config["install_dir"] + "biobakery_scripts/",
        main_dir=TMPDIR + HMN_OUTPUTS + MAIN
    shell:
        """
        {params.script_loc}get_counts_from_humann_logs.py --input {params.main_dir} --output {output}
        """

rule humann_feature_counts:
    input:
        genefam=rules.humann_join_relev.output.relev_genefam,
        pathabd=rules.humann_join_relev.output.relev_pathabd,
        regrouped=rules.humann_regroup.output 
    output:
        TMPDIR + HMN_OUTPUTS + COUNTS + "humann_feature_counts.tsv"
    params:
        script_loc=config["install_dir"] + "biobakery_scripts/",
        outdir=TMPDIR + HMN_OUTPUTS + COUNTS,
        regroups=config["name"]
    shell:
        """
        # counts of taxonomic features is currently not compatible with metaphlan 4
        # It isn't used for the humann_feature_counts, so I will remove it for now, but wanted to leave a note.

        {params.script_loc}count_features.py -i {input.genefam} -o {params.outdir}/genefamilies_counts.tsv --reduce-sample-name --ignore-un-features --ignore-stratification
        {params.script_loc}count_features.py -i {input.pathabd} -o {params.outdir}/pathabundance_counts.tsv --reduce-sample-name --ignore-un-features --ignore-stratification
        {params.script_loc}count_features.py -i {input.regrouped} -o {params.outdir}/{params.regroups}s_counts.tsv --reduce-sample-name --ignore-un-features --ignore-stratification

        humann_join_tables -i {params.outdir} -o {output} --file_name _counts.tsv
        """

rule metaphlan_biom:
    input:
        rules.combine_metaphlan.output.merged_file
    output:
        TMPDIR + MPHL_OUTPUTS + BIOMS + "metaphlan_taxonomic_profiles.biom"
    run:
        make_biom_from_tsv(input, output)

rule humann_pathabund_biom:
    input:
        rules.humann_join_relev.output.relev_pathabd
    output:
        TMPDIR + HMN_OUTPUTS + BIOMS + "humann_pathabundance_" + config["units"] + ".biom"
    run:
        make_biom_from_tsv(input, output)

rule humann_ecs_biom:
    input:
        rules.humann_rename.output
    output:
        TMPDIR + HMN_OUTPUTS + BIOMS + "humann_" + config["name"] + "_" + config["units"] + ".biom"
    run:
        make_biom_from_tsv(input, output)

rule visualizations:
    input:
        rules.humann_ecs_biom.output,
        rules.metaphlan_biom.output,
        rules.humann_pathabund_biom.output
    params:
        script_loc=config["install_dir"] + "biobakery_scripts/",
        yaml=os.environ.get("CONFIG_FILE"),
        file="visualizations.html",
        final_dir=TMPDIR + config["output_dir"],
        workdir=TMPDIR
    output:
        TMPDIR + config["output_dir"] + "visualizations.html"
    shell:
        """
        cp {params.script_loc}visualizations.qmd {params.final_dir}
        cp {params.yaml} {params.final_dir}

        cd {params.final_dir}

        quarto render visualizations.qmd --to html \
        -o {params.file} --execute-params {params.yaml}

        rm visualizations.qmd {params.yaml}
        cd {params.workdir}
        """
