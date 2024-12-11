rule analysis:
    input:
        "results/analysis/01_data-preparation.html",
        "results/analysis/02_normalization.html",
        "results/analysis/03_data-exploration.html",
        "results/analysis/04_expressed-genes.html",
        "results/analysis/05_functional-enrichment.html",
        "results/analysis/06_chromatin-states.html"
    default_target: True


rule data_preparation:
    input:
        qmd = "analysis/01_data-preparation.qmd",
        counts = "results/counts/all.tsv"
    output:
        report = "results/analysis/01_data-preparation.html",
        raw_counts = "results/R/01_data-preparation/raw_counts_SE.rds",
        coding_counts = "results/R/01_data-preparation/coding_counts_SE.rds"
    log:
        run = "logs/R/01_data-preparation.log",
        session = "logs/R/01_data-preparation_session.rds"
    conda: "P02-r"
    container: None
    shell:
        fmt("""
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """)


rule normalization:
    input:
        qmd = "analysis/02_normalization.qmd",
        counts = rules.data_preparation.output.coding_counts
    output:
        report = "results/analysis/02_normalization.html",
        normalized_counts = "results/R/02_normalization/normalized_counts_SE.rds"
    log:
        run = "logs/R/02_normalization.log",
        session = "logs/R/02_normalization_session.rds"
    conda: "P02-r"
    container: None
    shell:
        fmt("""
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """)


rule data_exploration:
    input:
        qmd = "analysis/03_data-exploration.qmd",
        normalized_counts = rules.normalization.output.normalized_counts
    output:
        report = "results/analysis/03_data-exploration.html"
    log:
        run = "logs/R/03_data-exploration.log",
        session = "logs/R/03_data-exploration_session.rds"
    conda: "P02-r"
    container: None
    shell:
        fmt("""
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """)


rule expressed_genes:
    input:
        qmd = "analysis/04_expressed-genes.qmd",
        normalized_counts = rules.normalization.output.normalized_counts
    output:
        report = "results/analysis/04_expressed-genes.html",
        expressions = "results/R/04_expressed-genes/expressions_cell_lines.rds",
        genesets = expand(
            "results/chromatin-states/annotation/{geneset}.txt",
            geneset = ["expressed-genes", "not-expressed-genes", "ignored-genes"]
        )
    log:
        run = "logs/R/04_expressed-genes.log",
        session = "logs/R/04_expressed-genes_session.rds"
    conda: "P02-r"
    container: None
    shell:
        fmt("""
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """)


rule functional_enrichment:
    input:
        qmd = "analysis/05_functional-enrichment.qmd",
        expressions = rules.expressed_genes.output.expressions
    output:
        report = "results/analysis/05_functional-enrichment.html"
    log:
        run = "logs/R/05_functional-enrichment.log",
        session = "logs/R/05_functional-enrichment.rds"
    conda: "P02-r"
    container: None
    shell:
        """
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """


rule chromatin_states_analysis:
    input:
        qmd = "analysis/06_chromatin-states.qmd",
        chromstates_enrichments = rules.run_chromatin_state_enrichments.input
    output:
        report = "results/analysis/06_chromatin-states.html"
    log:
        run = "logs/R/06_chromatin-states.log",
        session = "logs/R/06_chromatin-states_session.rds"
    conda: "P02-r"
    container: None
    shell:
        """
        quarto render {input.qmd} --to html > {log.run} 2>&1
        """
