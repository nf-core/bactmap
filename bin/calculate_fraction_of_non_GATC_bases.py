#!/usr/bin/env python3
from Bio import SeqIO
import argparse, sys, os

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print('{0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)

    def is_valid_file(self, parser, arg):
        if not os.path.isfile(arg):
            parser.error("The file %s does not exist!" % arg)
        else:
            return arg


def argparser():
    description = """
    A script to find the fraction of non-GATC bases in a fasta file
    """
    parser = ParserWithErrors(description = description)
    parser.add_argument("-f", "--fasta_file", required=True, 
                        help="fasta file path",
                        type=lambda x: parser.is_valid_file(parser, x))


    return parser


def calculate_fraction_of_non_GATC_bases(fasta_file):
    record  = SeqIO.read(fasta_file, 'fasta')
    total_len = len(record.seq)
    num_Gs = record.seq.upper().count('G')
    num_As = record.seq.upper().count('A')
    num_Ts = record.seq.upper().count('T')
    num_Cs = record.seq.upper().count('C')
    return (total_len - (num_Gs + num_As + num_Ts + num_Cs))/total_len

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()
    print(calculate_fraction_of_non_GATC_bases(args.fasta_file))
