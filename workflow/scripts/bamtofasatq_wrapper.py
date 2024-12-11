"""

Usage:
------
    input positional arguments:
        BAM file. In the case of paired-end reads the bam files need to be
        sorted first.
    output positional arguments:
        One (single-end reads) or Two (paired-end reads) .fastq.gz files in
        which the reads from the BAM file is going to be written.

Author: Markus Riedl
E-Mail: markus.riedl@boku.ac.at, markusriedl@acib.at
Created: 2023-01-31
"""


from os import remove
from os.path import basename, splitext
from typing import List
from snakemake.shell import shell
from pdb import set_trace


def is_paired(output: List[str]) -> bool:
    return len(output) > 1


def output_subshell(output: List[str], threads: int) -> str:
    if is_paired(output):
        outputcmd = (
            f"-fq  >(pigz -p {int(threads / 2)} -5 --stdout > {output[0]}) "
            f"-fq2 >(pigz -p {int(threads / 2)} -5 --stdout > {output[1]})"
        )
    else:
        outputcmd = f"-fq  >(pigz -p {threads} -5 --stdout > {output})"

    return outputcmd


outputcmd = output_subshell(snakemake.output, snakemake.threads)

if is_paired(snakemake.output):
    shell(
        "samtools sort"
        " -@ {snakemake.threads}"
        " -n"
        " -o input_sorted.bam"
        " {snakemake.input}"
        " |& tee -a {snakemake.log.run}"
    )
    bam_input = "input_sorted.bam"
else:
    bam_input = snakemake.input


shell(
    "bamToFastq"
    " -i {bam_input}"
    " {outputcmd}"
    " |& tee -a {snakemake.log.run}"
)


shell("""
cat <<-END_VERSION > {snakemake.log.version}
{snakemake.rule}:
  bedtools: "$(bedtools -version | grep -Eo [0-9.]+)"
  samtools: "$(samtools version | head -n1 | grep -Eo [0-9.]+)"
  date: "$(date -Iseconds)"
END_VERSION
"""
)
