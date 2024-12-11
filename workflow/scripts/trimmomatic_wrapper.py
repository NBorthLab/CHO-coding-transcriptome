"""
Custom wrapper for trimming of RNA-Seq samples in a snakemake rule using
trimmomatic.

Usage:
------
    input named arguments:
        reads:      (ONLY SINGLE-END) Raw, untrimmed single-end reads as
                    gzipped fastq.
        fw:         (ONLY PAIRED-END) Raw, untrimmed paired-end forward reads
                    as gzipped fastq.
        rv:         (ONLY PAIRED-END) Raw, untrimmed paired-end reverse reads
                    as gzipped fastq.
    output positional arguments:
        If single-end there is only one output as .fastq.gz.
        If paired-end there are four outputs (fw_paired, fw_unpaired,
            rv_paired, rv_unpaired).
    parameter argument:
        adapters:   Should be defined in the config yaml file and associate the
                    dataset to the corresponding adapters fasta.

Author: Markus Riedl
E-Mail: markus.riedl@boku.ac.at, markusriedl@acib.at
Created: 2023-01-27
"""


from os import remove
from os.path import basename, splitext
from typing import List
from snakemake.shell import shell
from snakemake.utils import Namedlist
from pdb import set_trace


def is_paired(input: Namedlist, output: Namedlist) -> bool:
    return len(input) > 1 and len(output) > 1


def decompress_input(*input: str) -> List[str]:
    """
    Returns a string to the decompressed fastq files.
    """
    tmp_files = [splitext(basename(file))[0] for file in input]
    for file, tmp_file in zip(input, tmp_files):
        shell(
            f"pigz -p {snakemake.threads} --decompress --stdout {file} "
            f"> {tmp_file}"
        )
    return tmp_files


def compress_output(tmp_files: List[str], output: List[str]) -> None:
    """
    Compresses the fastq file.
    """
    for file, tmp_file in zip(output, tmp_files):
        shell(
            f"pigz -p {snakemake.threads} -5 --stdout {tmp_file} > {file}"
        )


def cleanup(*files: str) -> None:
    """
    Clean up temporary files.
    """
    for file in files:
        print("Removing", file)
        remove(file)


# Execute ----------------------------------------

# Decompress .fastq.gz to .fastq
raw_fastqs: List[str] = decompress_input(*snakemake.input)

trimmed_fastqs: List[str] = [
    splitext(basename(file))[0] for file in snakemake.output
]

# Run trimming
end_parameter = "PE" if is_paired(snakemake.input, snakemake.output) else "SE"

adapter = "resources/adapters/" + snakemake.params.adapters

shell(
    "trimmomatic {end_parameter} "
    "-threads {snakemake.threads} "
    "-Xmx20000M "
    "{raw_fastqs} "
    "{trimmed_fastqs} "
    "ILLUMINACLIP:{adapter}:2:30:10 "
    "LEADING:15 "
    "TRAILING:15 "
    "SLIDINGWINDOW:4:15 "
    "MINLEN:36 "
    "|& tee {snakemake.log}"
)

compress_output(trimmed_fastqs, snakemake.output)

# set_trace()

shell("""
cat <<-END_VERSION > {snakemake.log.version}
{snakemake.rule}:
  trimmomatic: $(trimmomatic -version)
  date: "$(date -Iseconds)"
END_VERSION
""")

