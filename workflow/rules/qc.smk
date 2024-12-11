# vim: set syntax=snakemake:

from os.path import exists


# ================================================
#   FASTQC
# ================================================


def get_fastq(wildcards):
    if wildcards.step == "raw":
        ext = "fastq.gz"
    else:
        ext = "fq.gz"

    path = "results/reads/{step}/{dataset}/{file}.{ext}".format(
        step = wildcards.step,
        dataset = wildcards.dataset,
        file = wildcards.file,
        ext = ext
    )
    return path


rule fastqc:
    """ FastQC on the raw reads and on trimmed reads. """
    input:
        get_fastq
    output:
        "results/qc/{step}/{dataset}/{file}_fastqc.html",
        "results/qc/{step}/{dataset}/{file}_fastqc.zip"
    params:
        outdir = lambda wildcards, output: dirname(output[0])
    threads: 4
    log:
        run = "logs/fastqc/{step}/{dataset}/{file}.log",
        version = "logs/fastqc/{step}/{dataset}/{file}_versions.yaml"
    conda:
        "../envs/qc.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        fastqc --outdir {params.outdir} \\
                --threads {threads} \\
                --dir . \\
                {input} \\
            > {log.run} 2>&1

        cat <<-END_VERSION > {log.version}
        {rule}:
          fastqc: "$(fastqc --version | grep -oE [0-9.]+)"
          date: "$(date -Iseconds)"
          condaenv: "$CONDA_PREFIX"
        END_VERSION
        """)



# ================================================
#   MULTIQC
# ================================================


def get_dataset_runs(wildcards):
    """
    Input function for MultiQC. Includes all logs, reports etc. that multiqc
    can read for each step.
    That includes in raw reads the fastqc reports , and in trimmed reads the
    fastqc reports and salmon quant logs.
    """

    # Trimmed reads have a different naming scheme for paired-end fastq files.
    if wildcards.step == "raw":
        pairs = ["1", "2"]
    else:
        pairs = ["1_paired", "2_paired"]

    # Get run accession numbers of the dataset
    run_accessions = runs \
            .query(f"DatasetReadable == '{wildcards.dataset}'") \
            .RunAccession \
            .to_list()

    # Single-end and Paired-end reads have different file names.
    if is_paired_end(wildcards.dataset):
        file = expand(
            "{accession}_{pair}_fastqc",
            accession = run_accessions,
            pair = pairs
        )
    else:
        file = expand(
            "{accession}_fastqc",
            accession = run_accessions
        )

    paths = expand(
        "results/qc/{step}/{dataset}/{file}.{ext}",
        step = wildcards.step,
        dataset = wildcards.dataset,
        file = file,
        ext = ["html", "zip"]
    )

    # Include STAR alignment logs, if they exist...
    quant_logs = f"results/alignment/{wildcards.dataset}"
    if wildcards.step == "trimmed" and exists(quant_logs):
        paths.append(quant_logs)

    # Include featureCounts quantification summary, if it exists...
    count_summary = f"results/counts/{wildcards.dataset}"
    if wildcards.step == "trimmed" and exists(count_summary):
        paths.append(count_summary)

    return paths


rule multiqc:
    """ MultiQC for raw data before trimming """
    input:
        get_dataset_runs
    output:
        "results/qc/{step}/{dataset}/multiqc_report.html",
        directory("results/qc/{step}/{dataset}/multiqc_data")
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0])
    log:
        run = "logs/multiqc/{dataset}_{step}.log",
        version = "logs/multiqc/{dataset}_{step}_versions.yaml"
    conda:
        "../envs/qc.yaml"
    shell:
        fmt("""
        multiqc \\
                -o {params.outdir} \\
                {input} \\
            > {log.run} 2>&1

        cat <<-END_VERSION > {log.version}
        {rule}:
          multiqc: "$(multiqc --version | grep -oE [0-9.]+)"
          date: "$(date -Iseconds)"
          condaenv: "$CONDA_PREFIX"
        END_VERSION
        """)


rule test_multiqc:
    input:
        expand("results/qc/trimmed/{dataset}/multiqc_report.html",
               dataset = [
                   "ERP122753",
                   "SRP066848",
                   "SRP069883",
                   "SRP111368",
                   "SRP159459",
                   "SRP234382",
                   "SRP246348",
                   "SRP324587",
                   "INH01_Papez",
                   "INH02_NovakProducers",
                   "INH03_WeingunySubclones",
                   "INH04_WeingunySubclonability",
                   "INH05_Novak2019Clones",
                   "INH06_Ruckerbauer",
                   "INH07_NovakTempShift"
            ]
        )


rule raw_multiqc:
    input:
        expand(
            "results/qc/raw/{dataset}/multiqc_report.html",
            dataset = runs.DatasetReadable
        )

rule trimmed_multiqc:
    input:
        expand(
            "results/qc/trimmed/{dataset}/multiqc_report.html",
            dataset = runs.DatasetReadable.unique().tolist()
        )

rule run_multiqc:
    input:
        rules.raw_multiqc.input,
        rules.trimmed_multiqc.input
