# vim: set syntax=snakemake:


rule unzip_fasta:
    """
    --singularity-args "--bind /data/borth/mriedl/raw_data"
    """
    input:
        "resources/raw_data/genome/{genome}.fna.gz"
    output:
        "results/genome/{genome}.fna"
    threads: workflow.cores
    shell:
        "gunzip -c {input} > {output}"

use rule unzip_fasta as unzip_gtf with:
    input:
        "resources/raw_data/genome/{genome}.gtf.gz"
    output:
        "results/genome/{genome}.gtf"

use rule unzip_fasta as unzip_gff with:
    input:
        "resources/raw_data/genome/{genome}.gff.gz"
    output:
        "results/genome/{genome}.gff"
