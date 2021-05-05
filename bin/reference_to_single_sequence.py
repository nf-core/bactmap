#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
    A script to parse a filtered VCF and 
    """
    parser = ParserWithErrors(description = description)
    parser.add_argument("-r", "--reference_file", required=True, 
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-o", "--output_fasta_file", required=True,
                    help="file path to output fasta file")

    return parser

def combine_sequences(reference_sequence):
    records = list(SeqIO.parse(reference_sequence, "fasta"))
    new_sequence = ''.join([str(record.seq) for record in records])
    new_id = '-'.join([record.id for record in records])
    new_record = SeqRecord(Seq(new_sequence), id = new_id, description = '')
    return(new_record)


def write_sequence(filepath, record):
    with open(filepath, 'w') as output:
        SeqIO.write(record, output, "fasta")

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

    new_record = combine_sequences(args.reference_file)
    write_sequence(args.output_fasta_file, new_record)
