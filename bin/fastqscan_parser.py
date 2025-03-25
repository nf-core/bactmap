#!/usr/bin/env python

import pandas as pd
import json
import glob
import sys
import os
import argparse

def parser_args(args=None):
    """ 
    Function for input arguments for fastqscan_parser.py
    """
    Description = 'Collect fastq-scan outputs and create a summary table'
    Epilog = """Example usage: python fastqscan_parser.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-of", "--output_file", type=str, default="fastq-scan_summary.tsv", help="fastq-scan summary file (default: 'fastq-scan_summary.tsv').")
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

def json_to_dataframe(json_files):
    """ 
    Function to take list of json files and create a summary table
    """
    json_names = [i.replace('.json', '') for i in json_files]
    json_names_df = pd.DataFrame(json_names)
    json_names_df.columns = ['sample']
    jsons_data = {}

    for index, file in enumerate(json_files):
        with open(file, 'r') as f:
            json_text = json.loads(f.read())
            qc = json_text['qc_stats']
            jsons_data[index] = qc

    jsons_data_df = pd.DataFrame.from_dict(jsons_data, orient = 'index')
    json_merged_df = json_names_df.join(jsons_data_df)

    return json_merged_df

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of fastq-scan json files
    json_files = sorted(glob.glob('*.json'))

    ## Create dataframe
    json_df = json_to_dataframe(json_files)

    ## Write output file
    json_df.to_csv(args.output_file, sep = '\t', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())