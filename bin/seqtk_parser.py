#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import pandas as pd


def parser_args(args=None):
    """
    Function for input arguments for seqtk_parser.py
    """
    Description = "Collect seqtk comp outputs and create a summary table"
    Epilog = """Example usage: python seqtk_parser.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-of",
        "--output_file",
        type=str,
        default="mapping_summary.tsv",
        help="seqtk comp summary file (default: 'mapping_summary.tsv').",
    )
    return parser.parse_args(args)


def make_dir(path):
    """
    Function for making a directory from a provided path
    """
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def seqtk_to_dataframe(file_list):
    """
    Function for creating a dataframe from a list of seqtk comp files
    """
    seqtk_file_list = [pd.read_csv(f, sep="\t", header=None) for f in file_list]
    seqtk_df = pd.concat(seqtk_file_list, ignore_index=True)
    seqtk_df = seqtk_df.iloc[:, 0:6]
    seqtk_df.columns = ["sample", "ref_length", "#A", "#C", "#G", "#T"]
    column_list = ["#A", "#C", "#G", "#T"]
    seqtk_df["mapped"] = seqtk_df[column_list].sum(axis=1)
    seqtk_df["%ref mapped"] = seqtk_df["mapped"] / seqtk_df["ref_length"] * 100

    return seqtk_df


def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of seqtk comp tsv outputs
    seqtk_files = sorted(glob.glob("*.tsv"))

    ## Create dataframe
    seqtk_df = seqtk_to_dataframe(seqtk_files)

    ## Write output file
    seqtk_df.to_csv(args.output_file, sep="\t", header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
