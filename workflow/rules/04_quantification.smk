# vim: set syntax=snakemake:


# ================================================
#   Read Strandedness
# ================================================

rule infer_strandedness:
    input:
        annotation = (f"resources/raw_data/genome/"
                      f"{config['assembly']}_genomic.gtf.gz"),
        alignment = ("results/alignment/{dataset}/{accession}/"
                     "Aligned.sortedByCoord.out.bam"),
    output:
        strandedness = "results/counts/{dataset}/{accession}_strandedness.txt"
    conda:
        "../envs/rseqc.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        gunzip -c {input.annotation} > annotation.gtf
        gtf2bed < annotation.gtf > annotation.bed
        infer_experiment.py -r annotation.bed -i {input.alignment} \
            > infer_experiment.out.txt

        awk '
            BEGIN {{ FS=": "; paired_end = "false" }}
            /PairEnd/ {{ paired_end = "true" }}
            /\\".?\\+\\+/ {{ fw=$2 }}
            /\\".?\\+\\-/ {{ rev=$2 }}
            END {{
                if (fw > 0.75)
                    print "1\\nForward"
                else if (rev > 0.75)
                    print "2\\nReverse"
                else
                    print "0\\nUnstranded"

                if (paired_end == "true")
                    print "Paired End"
                else
                    print "Single End"

                print fw,rev
            }}
        ' infer_experiment.out.txt > {output.strandedness}
        """)


# ================================================
#   featureCounts
# ================================================


rule quantify_reads:
    """\
    Quantify/Summarize Reads (or Fragments) from the alignments using
    featureCounts.
    """
    input:
        annotation = (f"resources/raw_data/genome/"
                      f"{config['assembly']}_genomic.gtf.gz"),
        reads = ("results/alignment/{dataset}/{accession}/"
                 "Aligned.sortedByCoord.out.bam"),
        strandedness = rules.infer_strandedness.output.strandedness,
        _genome_removed = "logs/star_align/.genome-index-removed"
    output:
        counts = "results/counts/{dataset}/{accession}_counts.txt"
    log:
        run = "logs/featureCounts/{dataset}/{accession}.log",
        version = "logs/featureCounts/{dataset}/{accession}_versions.yaml"
    threads: 8
    conda: "../envs/featurecounts.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        gunzip -c {input.annotation} > annotation.gtf
        strand=$(head -n1 {input.strandedness})
        p=$(if grep Paired {input.strandedness} > /dev/null; then echo "-p"; fi)

        featureCounts \\
            -T {threads} \\
            ${{p}} \\
            -s ${{strand}} \\
            -a annotation.gtf \\
            -o {output.counts} \\
            {input.reads} \\
        |& tee {log.run}

        cat <<-END_VERSION > {log.version}
        {rule}:
          subread: "$(featureCounts -v 2>&1 | grep -Eo [0-9.]+)"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


rule test_quant:
    input:
        "results/counts/INH01_Papez/INH010100_counts.txt"


rule run_quantification:
    input:
        expand(
            "results/counts/{dataset}/{accession}_counts.txt",
            zip,
            dataset = samples["DatasetReadable"],
            accession = samples["SampleAccession"]
        )
    output:
        "results/counts/all.tsv"
    message: "Quantification finished!"
    shell:
        fmt("""
        input=({input})
        tail -n +2 ${{input[1]}} | cut -f 1 > {output}

        # The paste command can't modify in-place. Hence, a temporary file is
        # used.
        for file in ${{input[@]}}; do
            paste {output} <(tail -n +2 ${{file}} | cut -f 7) \
                > /dev/shm/tempfile
            cat /dev/shm/tempfile > {output}
        done

        rm /dev/shm/tempfile

        # Prettify sample names in header
        sed -i \
            -e 's/\/Aligned\.sortedByCoord\.out\.bam//g' \
            -e 's/results\/alignment\///g' \
            {output}
        """)

