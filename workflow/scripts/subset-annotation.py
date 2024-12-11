"""
Subset a given GTF annotation to genes specified in a gene list.
Output is the subset GTF as well as BED files containing individual features of
the genes:
    - gene
    - exons
    - TSS
    - TES
    - TSS +/- 2kb

Author: Markus Riedl
E-Mail: markus.riedl@boku.ac.at, markusriedl@acib.at
Created: November 2023
"""


from pygtftk.gtf_interface import GTF
from pathlib import Path, PurePath
from argparse import ArgumentParser
import logging as log

# %%

def parse_arguments():
    parser = ArgumentParser(
        prog = "Subset Annotation",
        epilog = (
            "The script produces a subset GTF annotation file containing only "
            "the genes specified in --genelist.\n"
            "Additionally, it exports individual features TSS, TES, genes and "
            "exons as well as a +/- 2kb window around the TSS into BED files."
        )
    )

    parser.add_argument("--annotation", help = "GTF annotation file")
    parser.add_argument("--genelist", help = "List of genes to subset")
    parser.add_argument("--log", help = "Log file",
                        default = "subset-annotation.log")
    parser.add_argument("--verbose", action = "store_true",
                        help = "Print log on screen")

    args = parser.parse_args()
    return args

# %%

def load_annotation(annotation_file):
    return GTF(annotation_file)

# %%

def subset_annotation(genelist_file, annotation):
    """\
    Subset a given GTF annotation to only include the genes specified in
    `genelist_file`. Individual features of the genes are exported to BED files:
        - coordinates of the full gene
        - coordinates of the exons
        - coordinate of the TSS (transcription start site)
        - coordinate of the TES (transcription end site)
        - coordinates of a 2kb window around the TSS (TSS +/- 2kb)
    """

    # Read annotation, subset by genelist
    with open(genelist_file) as genelist:
        gl_annot = annotation.select_by_key(
            key              = "gene_id",
            file_with_values = genelist
        )

    # Output file prefix
    genelist_file = Path(genelist_file)
    fileprefix = str(genelist_file.parent) + "/" + genelist_file.stem

    # Export GTF annotation of the gene list
    with open(fileprefix + ".gtf", "w") as f:
        gl_annot.write(f)
        log.info("Saving subset annotation for genelist %s as GTF to %s",
                 genelist_file, fileprefix + ".gtf")

    # Get individual feature annotations for the genelist
    gl_annot_genes = gl_annot.select_by_key("feature", "gene")
    gl_annot_exons = gl_annot.select_by_key("feature", "exon")
    gl_annot_tss = gl_annot.get_tss()
    gl_annot_tes = gl_annot.get_tts()

    # Save feature annotations to BED
    with open(fileprefix + "_genes.bed", "w") as f:
        gl_annot_genes.write_bed(f)

    with open(fileprefix + "_exons.bed", "w") as f:
        gl_annot_exons.write_bed(f)

    gl_annot_tss.saveas(fileprefix + "_tss.bed")

    gl_annot_tes.saveas(fileprefix + "_tes.bed")

    # TSS +- 2kb and saving to BED
    with open(fileprefix + "_tss2kb.bed", "w") as file:

        for tss in gl_annot_tss:
            tss_upstream = tss.start - 2_000
            tss_downstream = tss.end + 2_000

            feature_line = [
                tss.chrom,
                str(tss_upstream),
                str(tss_downstream),
                tss.name,
                tss.score,
                tss.strand
            ]

            file.write("\t".join(feature_line) + "\n")

    log.info("Exported all feature BED files. Good bye")

# %%

def main():
    args = parse_arguments()

    log_handlers = [log.FileHandler(args.log, mode = "w")]
    if args.verbose:
        log_handlers.append(log.StreamHandler())

    log.basicConfig(
        level = log.INFO,
        force = True,
        format = "%(asctime)s :: %(levelname)s - %(message)s",
        datefmt = "%Y-%m-%d %H:%M:%S%z",
        handlers = log_handlers
    )

    genelist_file = args.genelist
    annotation_file = args.annotation

    log.info((
        f"Subsetting annotataion with the following parameters:\n"
        f"    genelist: {genelist_file}\n"
        f"    annotation: {annotation_file}"
        f"    log: {args.log}"
    ))

    annotation = load_annotation(annotation_file)
    subset_annotation(genelist_file, annotation)


if __name__ == "__main__":
    main()
