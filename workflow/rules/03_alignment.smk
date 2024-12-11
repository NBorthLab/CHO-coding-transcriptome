# vim: set syntax=snakemake:


rule star_index:
    input:
        fasta = expand(
            "resources/raw_data/genome/{assembly}_genomic.fna.gz",
            assembly = config["assembly"]
        ),
        gtf = expand(
            "resources/raw_data/genome/{assembly}_genomic.gtf.gz",
            assembly = config["assembly"]
        )
    output:
        directory("results/genome_index")
    threads: workflow.cores
    conda:
        "../envs/star.yaml"
    log:
        run = "logs/star_index.log",
        version = "logs/star_index_versions.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        gunzip -d -c {input.fasta} > genome.fna &
        fasta_pid=$!
        gunzip -d -c {input.gtf} > genome.gtf &
        gtf_pid=$!
        wait $fasta_pid $gtf_pid

        STAR \\
                --runThreadN {threads} \\
                --runMode genomeGenerate \\
                --genomeFastaFiles genome.fna \\
                --sjdbGTFfile genome.gtf \\
                --outTmpDir /dev/shm/STARtmp \\
                --genomeDir {output} \\
            |& tee {log.run}

        cat <<-END_VERSION > {log.version}
        {rule}:
          star: "$(STAR --version)"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


def get_reads_to_align(wildcards):
    """\
    Returns the corresponding path depending on whether the reads had to be
    trimmed or not.
    """

    # Get the run accession numbers corresponding to the RNA-seq sample
    corresponding_runs = runs[runs.SampleAccession == wildcards.accession] \
            ["RunAccession"]

    if is_paired_end(wildcards.dataset):
        return {
            "reads1": expand(
                f"results/reads/trimmed/{wildcards.dataset}/{{run}}_1_paired.fq.gz",
                run = corresponding_runs
            ),
            "reads2": expand(
                f"results/reads/trimmed/{wildcards.dataset}/{{run}}_2_paired.fq.gz",
                run = corresponding_runs
            )
        }
    else:
        return {
            "reads1": expand(
                f"results/reads/trimmed/{wildcards.dataset}/{{run}}.fq.gz",
                run = corresponding_runs
            )
        }


rule star_align:
    """\
    Align reads using STAR and an index of the genome. Loads the genome index
    into memory only once, will be removed after the last job finshed.

    Requires setting the ulimit -n (open files) to max!
    i.e. ulimit -n $(ulimit -Hn)
    """
    input:
        unpack(get_reads_to_align),
        index = rules.star_index.output
    output:
        alignment = "results/alignment/{dataset}/{accession}/Aligned.sortedByCoord.out.bam",
        reads_per_gene = "results/alignment/{dataset}/{accession}/ReadsPerGene.out.tab",
        log = "results/alignment/{dataset}/{accession}/Log.out",
        log_progress = "results/alignment/{dataset}/{accession}/Log.progress.out",
        log_final = "results/alignment/{dataset}/{accession}/Log.final.out"
    threads: workflow.cores
    log:
        run = "logs/star_align/{dataset}/{accession}.log",
        version = "logs/star_align/{dataset}/{accession}_versions.yaml"
    shadow: "shallow"
    conda: "../envs/star.yaml"
    script:
        "../scripts/star_wrapper.py"


rule index_alignment:
    """\
    Create .bai index of alignment .bam files. This is needed for the JBrowse
    Genome Browser.
    """
    input:
        rules.star_align.output.alignment
    output:
        "results/alignment/{dataset}/{accession}/Aligned.sortedByCoord.out.bam.bai"
    log:
        run = "logs/star_align/{dataset}/{accession}_alignment_index.log",
        version = "logs/star_align/{dataset}/{accession}_alignment_index_versions.yaml"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        fmt("""
        samtools index \\
            -@ {threads} \\
            -o {output} \\
            {input} \\
        > {log.run} 2>&1

        cat <<-END_VERSION > {log.version}
        {rule}:
          star: "$(samtools version | head -n 1 | grep -oE [0-9.]+)"
          date: "$(date -Iseconds)"
          condaenv: "$CONDA_PREFIX"
        END_VERSION
        """)


rule test_alignment:
    input:
        "results/alignment/INH01_Papez/INH010100/Aligned.sortedByCoord.out.bam"


rule debug_alignment:
    input:
        "results/alignment/INH02_NovakProducers/INH020700/Aligned.sortedByCoord.out.bam"


rule run_alignment:
    input:
        expand(
            "results/alignment/{dataset}/{accession}/Aligned.sortedByCoord.out.bam",
            zip,
            dataset = samples["DatasetReadable"],
            accession = samples["SampleAccession"]
        ),
        index = rules.star_index.output
    output:
        touch("logs/star_align/.genome-index-removed")
    message: "STAR alignment finished. Removing genome index from memory."
    conda: "../envs/star.yaml"
    shell:
        """
        STAR --genomeDir {input.index} --genomeLoad Remove
        """


rule run_index_alignment:
    input:
        expand(
            "results/alignment/{dataset}/{accession}/Aligned.sortedByCoord.out.bam.bai",
            zip,
            dataset = samples["DatasetReadable"],
            accession = samples["SampleAccession"]
        )

