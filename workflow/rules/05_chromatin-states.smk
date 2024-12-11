import textwrap


rule create_annotations_files:
    """\
    Create annotation files for genes, TSS, TES, exons, TSS-2kb

    NOTE: Not a generalized rule. The PICRH_aliases file is specific to the
    PICRH genome assembly.
    """
    input:
        annotation = (
            f"resources/raw_data/genome/{config['assembly']}_genomic.gtf.gz"
        ),
        genelist = "results/chromatin-states/annotation/{genelist}.txt",
        aliases = "resources/PICRH_aliases.txt"
    output:
        multiext(
            "results/chromatin-states/annotation/{genelist}",
            ".gtf", "_genes.bed", "_exons.bed", "_tss.bed", "_tes.bed", "_tss2kb.bed"
        )
    log:
        run = "logs/chromatin-states/subset-annotation_{genelist}.log",
        version = "logs/chromatin-states/subset-annotation_{genelist}_versions.yaml"
    conda: "../envs/python.yaml"
    shadow: "shallow"
    shell:
        fmt("""
        gunzip -d -c {input.annotation} > annotation.gtf

        # Convert the chromosome identifiers in the GTF annotation to be
        # compatible with the chromatin states results.
        # First, swap the columns in the alias file.
        # That same alias file was used to rename the chromosomes for the
        # chromatin states estimation workflow in P01-GenomeBrowser (but the
        # other way around of course).
        awk '{{ print $2 "\\t" $1 }}' {input.aliases} > aliases.txt

        # Now, Substitute the names in the first column of the gtf files with
        # the format used in the chromatin states results.
        sed -f <(printf 's/%s/%s/g\\n' $(< aliases.txt)) annotation.gtf \\
                > annotation_aliased.gtf

        python workflow/scripts/subset-annotation.py \\
            --annotation annotation_aliased.gtf \\
            --genelist {input.genelist} \\
            --log {log.run}

        cat <<-END_VERSION > {log.version}
        {rule}:
          python: "$(python --version | grep -oE [0-9.]+)"
          pygtftk: "$(conda list | awk '/pygtftk/ {{ print $2 }}')"
          date: "$(date -Iseconds)"
          condaenv: "$CONDA_PREFIX"
        END_VERSION
        """)


rule run_annotation_subsetting:
    input:
        "results/chromatin-states/annotation/expressed-genes_tss.bed",
        "results/chromatin-states/annotation/not-expressed-genes_tss.bed",
        "results/chromatin-states/annotation/ignored-genes_tss.bed"


rule overlap_enrichment:
    """\
    Perform overlap enrichment analysis on Chromatin States with regions of
    expressed/not expressed genes (gene, exons, TSS +/- 2kb window).
    """
    input:
        genomic_regions = expand(
            "results/chromatin-states/annotation/{{expression}}_{region}.bed",
            region = ["genes", "exons", "tss2kb"]
        ),
        chromatin_states = "resources/raw_data/chromatin-states/{tp}"
    output:
        "results/chromatin-states/overlap-enrichment/{expression}_{tp}.txt"
    conda: "../envs/chromhmm.yml"
    container: None # ChromHMM requires X11 to create images
    log:
        run = "logs/chromatin-states/overlap-enrichment/{expression}_{tp}.log",
        version = "logs/chromatin-states/overlap-enrichment/{expression}_{tp}_versions.yaml"
    params:
        outprefix = lambda wildcards, output: splitext(output[0])[0],
        regionsdir = lambda wildcards, input: dirname(input.genomic_regions[0])
    threads: 2
    shadow: "shallow"
    shell:
        fmt("""
        REGIONS=({input.genomic_regions})
        # Get only the filenames, without the path.
        echo "${{REGIONS[@]##*/}}" | sed 's/ /\\n/g' > coordsfile.txt

        ChromHMM.sh -Xms2G -Xmx250G \\
                OverlapEnrichment \\
                -posterior \\
                -f coordsfile.txt \\
                {input.chromatin_states} \\
                {params.regionsdir} \\
                {params.outprefix} \\
            > {log.run} 2>&1

        cat <<-END_VERSION > {log.version}
        {rule}:
          chromhmm: "$(ChromHMM.sh Version | sed -E 's/.*Version.([0-9.]+).*/\\1/')"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


rule neighborhood_enrichment:
    """\
    Perform neighborhood enrichment analysis on Chromatin States with positions
    of expressed/not expressed genes (TSS and TES).
    """
    input:
        genomic_positions = "results/chromatin-states/annotation/{position}.bed",
        chromatin_states = "resources/raw_data/chromatin-states/{tp}"
    output:
        "results/chromatin-states/neighborhood-enrichment/{position}_{tp}.txt"
    conda: "../envs/chromhmm.yml"
    container: None # ChromHMM requires X11 to create images
    log:
        run = "logs/chromatin-states/neighborhood-enrichment/{position}_{tp}.log",
        version = "logs/chromatin-states/neighborhood-enrichment/{position}_{tp}_versions.yaml"
    params:
        outprefix = lambda wildcards, output: splitext(output[0])[0]
    threads: 2
    shadow: "shallow"
    shell:
        fmt("""
        ChromHMM.sh -Xms2G -Xmx250G \\
                NeighborhoodEnrichment \\
                -posterior \\
                -colfields 0,1,5 \\
                {input.chromatin_states} \\
                {input.genomic_positions} \\
                {params.outprefix} \\
            > {log.run} 2>&1

        cat <<-END_VERSION > {log.version}
        {rule}:
          chromhmm: "$(ChromHMM.sh Version | sed -E 's/.*Version.([0-9.]+).*/\\1/')"
          date: "$(date -Iseconds)"
        END_VERSION
        """)


rule run_chromatin_state_enrichments:
    input:
        expand(
            "results/chromatin-states/overlap-enrichment/{expression}_{tp}.txt",
            expression = [
                "expressed-genes",
                "not-expressed-genes",
                "ignored-genes",
                "reactive-genes",
                "non-reactive-genes"
            ],
            tp = ["Tp5", "Tp6", "Tp7", "Tp8", "Tp9", "Tp10", "Tp11", "Tp12"]
        ),
        expand(
            "results/chromatin-states/neighborhood-enrichment/{genelist}_{position}_{tp}.txt",
            genelist = [
                "expressed-genes",
                "not-expressed-genes",
                "ignored-genes",
                "reactive-genes",
                "non-reactive-genes"
            ],
            position = ["tss", "tes"],
            tp = ["Tp5", "Tp6", "Tp7", "Tp8", "Tp9", "Tp10", "Tp11", "Tp12"]
        )
