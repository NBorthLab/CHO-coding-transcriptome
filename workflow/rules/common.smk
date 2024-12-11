# vim: set ft=python:


from typing import Dict
import pandas as pd
import yaml
import re
from textwrap import dedent


runs: pd.DataFrame = pd.read_csv(config["runs_csv"])

# Subset runs according to Source
ncbi_runs: pd.DataFrame = runs.query("Source == 'SRA'")
inhouse_runs: pd.DataFrame = runs.query("Source == 'In-House'")

samples: pd.DataFrame = runs.drop(columns = "RunAccession").drop_duplicates()

datasets: pd.DataFrame = samples.drop(columns = "SampleAccession").drop_duplicates()


def is_paired_end(dataset: str) -> bool:
    """ Whether the dataset is Paired-end (True) or Single-end (False). """
    return datasets.query(f"DatasetReadable == '{dataset}'").PairedEnd.iloc[0]


# In-House CSV needed to access FilePaths of original BAM files
inhouse_data: pd.DataFrame = pd.read_csv("resources/inhouse_samples.csv") \
    .query("DatasetAccession != 'INH080000'")

test_set: pd.DataFrame = runs \
    .query(
        "SampleAccession in ['INH010100', 'ERS4805133']"
    )

inhouse_data.loc[inhouse_data.RunAccession == 'INH010101', "FilePath"].iloc[0]


def fmt(text):
    """
    Prettify shell commands for logs in Snakemake rules
    """
    return dedent(text)
