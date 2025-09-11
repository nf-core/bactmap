#!/usr/bin/env python
## Originally written by Andries van Tonder and released under the MIT license.
## See git repository (https://github.com/nf-core/bactmap) for full license text.

import pandas as pd
import os
import sys
import glob
import json
import argparse

def parser_args(args=None):
    """
    Function for input arguments for read_stats.py
    """
    Description = 'Collect fastq-scan and create a table for each sample'
    Epilog = """Example usage: python read_stats.py """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-of", "--output_file" , type=str, default="read_stats.csv", help="read stats file (default: 'read_stats.csv').")
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
    json_names_df.columns = ['Sample']
    jsons_data = {}

    for index, file in enumerate(json_files):
        with open(file, 'r') as f:
            json_text = json.loads(f.read())
            qc = json_text['qc_stats']
            jsons_data[index] = qc

    jsons_data_df = pd.DataFrame.from_dict(jsons_data, orient = 'index')
    json_merged_df = json_names_df.join(jsons_data_df)
    json_merged_df = json_merged_df.iloc[:, 0:4]

    return json_merged_df

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Create list of raw reads fastq-scan json files
    raw_json_files = sorted(glob.glob('*.raw.json'))

    ## Create list of processed reads fastq-scan json files
    processed_json_files = sorted(glob.glob('*.processed.json'))

    ## Create dataframe of raw reads fastq-scan results
    raw_json_df = json_to_dataframe(raw_json_files)
    raw_json_df = raw_json_df.rename(columns = {'total_bp' : 'raw_total_bp', 'read_total' : 'num_raw_reads', 'coverage' : 'raw_coverage'})
    raw_json_df['Sample'] = raw_json_df['Sample'].str.replace('.raw','')

    ## Create dataframe of processed reads fastq-scan results
    processed_json_df = json_to_dataframe(processed_json_files)
    processed_json_df = processed_json_df.rename(columns = {'total_bp' : 'processed_total_bp', 'read_total' : 'num_processed_reads', 'coverage' : 'processed_coverage'})
    processed_json_df['Sample'] = processed_json_df['Sample'].str.replace('.processed','')

    ## Merge fastq-scan dataframes
    fastqscan_merged = pd.merge(raw_json_df, processed_json_df, on = ['Sample'])
    fastqscan_merged['%reads_after_processed'] = fastqscan_merged['num_processed_reads'] / fastqscan_merged['num_raw_reads'] * 100

    ## Write output file
    fastqscan_merged.to_csv(args.output_file, sep = ',', header = True, index = False)

if __name__ == '__main__':
    sys.exit(main())
