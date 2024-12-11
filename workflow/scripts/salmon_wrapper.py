"""
Custom wrapper for salmon quant in mapping-based mode.

Usage:
------
    inputs:
        reads:      (SINGLE-END ONLY) Path to trimmed .fq.gz reads.
        r1, r2:     (PAIRED-END ONLY) Path to trimmed .fq.gz reads (forward and
                    reverse reads, respectively).
        index:      Paths to the salmon index files.
    output:
        quant:      Path for quant.sf output file.

Author: Markus Riedl
E-Mail: markus.riedl@boku.ac.at, markusriedl@acib.at
Created: 2023-02-06
"""

from os.path import basename, splitext, dirname
from typing import List
from snakemake.shell import shell
from pdb import set_trace


# Paired-end samples are expected to have 3 inputs, Single-end only two.
paired = "r1" in snakemake.input.keys()

index_directory = dirname(snakemake.input.index[0])

if paired:
    input_line = f"-1 {snakemake.input.r1} -2 {snakemake.input.r2}"
else:
    input_line = f"-r {snakemake.input.reads}"

output_directory = dirname(snakemake.output.quant)

shell(
    "salmon quant "
        "--threads {snakemake.threads} "
        "-i {index_directory} "
        "-l A "
        "{input_line} "
        "--validateMappings "
        "-o {output_directory} "
    "|& tee {snakemake.log}"
)
