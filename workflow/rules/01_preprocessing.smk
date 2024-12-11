# vim: set syntax=snakemake:


# =============================================================================
#   Published datasets from NCBI
# =============================================================================


# ================================================
#   SRA CONVERSION TO FASTQ
# ================================================


# Prefetch SRA -----------------------------------

rule fetch_sra:
    """\
    Fetch SRA from NCBI.

    NOTE: Should be made a temporary output later.

    For explanation of the regex see https://regex101.com/r/wZxS6b/1
    """
    output:
        directory("resources/raw_data/sra/{accession}")
    log: "logs/prefetch/{accession}_versions.yaml"
    retries: 3
    conda:
        "../envs/tools.yaml"
    shell:
        fmt("""
        prefetch {wildcards.accession} -O {output}

        cat <<-END_VERSION > {log}
        {rule}:
          sra_tools: "$(prefetch --version | grep -Eo [0-9.]+)"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


# Convert to .fastq -------------------------------


rule sra2fastq_PE:
    """\
    Convert from NCBI SRA to .fastq.gz.
    """
    input:
        rules.fetch_sra.output
    output:
        "results/reads/raw/{dataset}/{accession}_1.fastq.gz",
        "results/reads/raw/{dataset}/{accession}_2.fastq.gz"
    params:
        outdir = lambda wildcards, output: dirname(output[0])
    wildcard_constraints:
        dataset = "SRP.*|ERP.*",
        accession = "[A-Z0-9]+"
    log:
        run = "logs/sra2fastq/{dataset}_{accession}.log",
        version = "logs/sra2fastq/{dataset}_{accession}_versions.yaml"
    threads: 6
    retries: 3
    conda:
        "../envs/tools.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        fasterq-dump \\
                -e {threads} \\
                --force \\
                --split-3 \\
                ./{input}/{wildcards.accession} \\
            |& tee {log.run}

        pigz -p {threads} *.fastq
        mv *.fastq.gz {params.outdir}

        cat <<-END_VERSION > {log.version}
        {rule}:
          sra_tools: "$(fasterq-dump --version | grep -Eo [0-9.]+)"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


use rule sra2fastq_PE as sra2fastq_SE with:
    output:
        "results/reads/raw/{dataset}/{accession}.fastq.gz",


rule run_sra2fastq:
    """\
    Runs the the conversion from the binary SRA to gzipped FASTQ files.
    Separates the datasets that require merging before further processing.
    """
    input:
        expand(
            "results/reads/raw/{dataset}/{accession}.fastq.gz",
            zip,
            dataset = ncbi_runs["DatasetAccession"],
            accession = ncbi_runs["RunAccession"]
        )


# =============================================================================
#   In-House datasets
# =============================================================================


# ================================================
#   BAM to FASTQ
# ================================================

# Symlink BAM files ------------------------------

def get_bam_filepaths(wildcards):
    """ Return the server file path of the in-house run accession. """
    # accession = inhouse_runs["RunAccession"] == wildcards.accession
    return inhouse_data \
        .loc[inhouse_data.RunAccession == wildcards.accession, "FilePath"] \
        .iloc[0]


rule symlink_bam:
    """\
    Create a symbolic link to the raw data from the various directories on the
    server.
    """
    output:
        "resources/raw_data/bam/{datasets}/{accession}.bam"
    params:
        source_file = get_bam_filepaths
    container: None
    shell:
        """
        ln -s {params.source_file} {output}
        """


rule test_symlink_bam:
    input:
        "resources/raw_data/bam/INH01_Papez/INH010101.bam"


rule run_symlink_bam:
    input:
        expand("resources/raw_data/bam/{datasets}/{accession}.bam",
            zip,
            datasets = inhouse_runs["DatasetReadable"],
            accession = inhouse_runs["RunAccession"]
        )


# Sort BAM files ---------------------------------

rule sort_bam:
    """\
    --- DEPRECATED ---

    Sort BAM files of paired-end reads.
    """
    input:
        "resources/raw_data/bam/{dataset}/{accession}.bam"
    output:
        temp("results/reads/bam/{dataset}/{accession}_sort.bam")
    log:
        run = "logs/sort_bam/{dataset}_{accession}.log",
        version = "logs/sort_bam/{dataset}_{accession}_versions.yaml"
    conda: "../envs/samtools.yaml"
    threads: 4
    shell:
        """
        samtools sort \\
                -@ {threads} \\
                -T /dev/shm \\
                -n \\
                -o {output} \\
                {input} \\
            |& tee {log.run}

        cat <<-END_VERSION > {log.version}
        {rule}:
          samtools: "$(samtools version | head -n1 | grep -Eo [0-9.]+)"
          date: "$(date -Iseconds)"
        END_VERSION
        """


# Convert to .fastq ------------------------------

rule bamtofastq_SE:
    """\
    Convert BAM to fastq files.

    Dataset regex explanation: https://regex101.com/r/UtkejB/1
    """
    input:
        "resources/raw_data/bam/{dataset}/{accession}.bam"
    output:
        "results/reads/raw/{dataset}/{accession}.fastq.gz"
    wildcard_constraints:
        accession = "[A-Z0-9]+",
        dataset = "INH.*"
    log:
        run = "logs/bamtofastq/{dataset}_{accession}.log",
        version = "logs/bamtofastq/{dataset}_{accession}_versions.yaml"
    threads: 8
    shadow: "shallow"
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/bamtofasatq_wrapper.py"


use rule bamtofastq_SE as bamtofastq_PE with:
    output:
        "results/reads/raw/{dataset}/{accession}_1.fastq.gz",
        "results/reads/raw/{dataset}/{accession}_2.fastq.gz"


rule run_bamtofastq:
    input:
        expand(
            "results/reads/raw/{dataset}/{accession}.fastq.gz",
            zip,
            dataset = inhouse_runs.loc[~inhouse_runs.PairedEnd, "DatasetReadable"],
            accession = inhouse_runs.loc[~inhouse_runs.PairedEnd, "RunAccession"]
        ),
        expand(
            expand(
                "results/reads/raw/{dataset}/{accession}_{{num}}.fastq.gz",
                zip,
                dataset = inhouse_runs.loc[inhouse_runs.PairedEnd, "DatasetReadable"],
                accession = inhouse_runs.loc[inhouse_runs.PairedEnd, "RunAccession"]
            ),
            num = ["1", "2"]
        )

