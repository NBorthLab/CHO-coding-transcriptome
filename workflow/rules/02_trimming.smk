# vim: set syntax=snakemake:


# ================================================
#   TRIMMING
# ================================================


rule trim_samples_PE:
    """ Some samples require trimming. """
    input:
        fw = "results/reads/raw/{dataset}/{accession}_1.fastq.gz",
        rv = "results/reads/raw/{dataset}/{accession}_2.fastq.gz"
    output:
        "results/reads/trimmed/{dataset}/{accession}_1_paired.fq.gz",
        "results/reads/trimmed/{dataset}/{accession}_1_unpaired.fq.gz",
        "results/reads/trimmed/{dataset}/{accession}_2_paired.fq.gz",
        "results/reads/trimmed/{dataset}/{accession}_2_unpaired.fq.gz"
    threads: 8
    log:
        run = "logs/trimmomatic/{dataset}_{accession}.log",
        version = "logs/trimmomatic/{dataset}_{accession}_versions.yaml"
    params:
        adapters = lambda wildcards: config["adapters"][wildcards.dataset]
    shadow: "shallow"
    conda:
        "../envs/trimmomatic.yaml"
    script:
        "../scripts/trimmomatic_wrapper.py"


use rule trim_samples_PE as trim_samples_SE with:
    input:
        reads = "results/reads/raw/{dataset}/{accession}.fastq.gz"
    output:
        "results/reads/trimmed/{dataset}/{accession}.fq.gz"
    threads: 8


rule test_wrapper_trim:
    input:
        "results/reads/raw/INH03_WeingunySubclones/INH030101.fastq.gz"
    output:
        "tmp/trim_test.fastq.gz"
    log: "tmp/trim_test.log"
    params:
        trimmer = ["TRAILING:3"]
    threads: 20
    resources:
        mem_mb = 20000
    wrapper:
        "v1.21.4/bio/trimmomatic/se"


rule test_trim:
    input:
        "results/reads/trimmed/INH01_Papez/INH010101_1_paired.fq.gz"


rule run_trimming:
    input:
        expand(
            "results/reads/trimmed/{dataset}/{accession}.fq.gz",
            zip,
            dataset = runs[~runs.PairedEnd]["DatasetReadable"],
            accession = runs[~runs.PairedEnd]["RunAccession"]
        ),
        expand(
            "results/reads/trimmed/{dataset}/{accession}_1_paired.fq.gz",
            zip,
            dataset = runs[runs.PairedEnd]["DatasetReadable"],
            accession = runs[runs.PairedEnd]["RunAccession"]
        )
