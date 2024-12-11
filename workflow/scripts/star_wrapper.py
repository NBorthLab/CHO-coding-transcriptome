"""
Custom wrapper for STAR alignment.

Input
-----
reads1 : list of files
    List of (forward) reads to align to the genome. Paired- or Single-end.
    Might be length of one in case there are no technical replicates.
reads2 : list of files, optional
    List of reverse reads to align to the genome. Only with paired-end reads.
    Must be same length as reads1 and might be just one (depending on technical
    replicates).
index : directory
    Path to the directory containing the genome index.
annotation : file
    Path to the .gtf or .gff annotation file.

Output
------
alignment : file
    Sorted BAM file of the alignment.
genecounts : file
    Gene counts from STAR.
log, log_final, log_progress : file
    Some logs.


Author: Markus Riedl
E-Mail: markus.riedl@boku.ac.at, markusriedl@acib.at
Created: 2023-02-08
"""


from os import remove
from os.path import basename, splitext
from tempfile import TemporaryDirectory
from typing import List, Dict
from snakemake.shell import shell
from snakemake.utils import Namedlist
from pdb import set_trace


# Type definitions ------------------------------

InputDict = Dict[str, List[str]]


# Functions --------------------------------------

def format_inputs(input: InputDict) -> str:
    """
    Format command line string for inputs to STAR. Considers whether the data
    is paired- or single-end and also if there are technical replicates or not.
    """
    if "reads2" in input.keys():
        return ",".join(input["reads1"]) + " " + ",".join(input["reads2"])
    else:
        return ",".join(input["reads1"])



# Execution --------------------------------------

input_files = format_inputs(snakemake.input)


shell(
    "STAR"
    " --runThreadN {snakemake.threads}"
    " --genomeDir {snakemake.input.index}"
    " --readFilesCommand gunzip -c"
    " --readFilesIn {input_files}"
    " --genomeLoad LoadAndKeep"
    " --quantMode GeneCounts"
    " --outSAMtype BAM SortedByCoordinate"
    " --limitBAMsortRAM 256000000000"
    " --outTmpDir STARtmp"
    " --outFileNamePrefix ./"
    " |& tee {snakemake.log}"
)

shell("""
    copy_jobs=()
    cp ReadsPerGene.out.tab {snakemake.output.reads_per_gene} &
    copy_jobs+=($!)
    cp Log.final.out {snakemake.output.log_final} &
    copy_jobs+=($!)
    cp Log.progress.out {snakemake.output.log_progress} &
    copy_jobs+=($!)
    cp Log.out {snakemake.output.log} &
    copy_jobs+=($!)
    cp Aligned.sortedByCoord.out.bam {snakemake.output.alignment} &
    copy_jobs+=($!)
    wait ${{copy_jobs[@]}}
""")

shell("""
cat <<-END_VERSION > {snakemake.log.version}
{snakemake.rule}:
  star: $(STAR --version)
  date: "$(date -Iseconds)"
END_VERSION
""")
