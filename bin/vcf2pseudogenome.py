#!/usr/bin/env python
from pysam import VariantFile
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
    parser.add_argument("-b", "--filtered_bcf_file", required=True,
                        help="filtered bcf file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-o", "--output_fasta_file", required=True,
                    help="file path to output fasta file")

    return parser

def calculate_reference_length(reference_file):
    record = SeqIO.read(reference_file, format = 'fasta')
    return len(record.seq)

def calculate_gaps_to_add(gap_start_position, gap_end_position):
    return '-' * (gap_end_position - gap_start_position)

def filtered_bcf_to_fasta(filtered_bcf_file, reference_length):
    with VariantFile(filtered_bcf_file) as vcf_reader:
        previous_pos = 0
        sequence = ""
        for record in vcf_reader.fetch():
            if record.pos == previous_pos: # Insertion - remove previous character and add 'N'
                sequence = sequence[:-1] # remove last position
                sequence += 'N'
            else:
                if previous_pos + 1 < record.pos: # large deletion add gaps before adding next position
                    sequence += calculate_gaps_to_add(previous_pos + 1, record.pos)
                if 'PASS' in record.filter.keys(): # HQ SNP
                    gt = record.samples[0]['GT'][0] # get genotype (1st from tuple)
                    if 'INDEL' in record.info: # indel
                        sequence += 'N'
                    elif gt == 0: # reference base
                        sequence += record.ref.lower() # add reference base as lower case
                    else: # variant
                        if len(record.alts) != 1: # if more than one ALT genotype so add N
                            sequence += 'N'
                        else: # add ALT SNP as upper case
                            sequence += record.alts[gt-1].upper()
                else: # if not PASS it's a low qual SNP so add N
                    sequence += 'N'
            previous_pos = record.pos
        if previous_pos != reference_length: # if gap at the end
            sequence += calculate_gaps_to_add(previous_pos, reference_length)
        return sequence

def write_sequence(filepath, id, sequence):
    with open(filepath, 'w') as output:
        record = SeqRecord(Seq(sequence), id = id, description = '')
        SeqIO.write(record, output, "fasta")

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

    reference_length = calculate_reference_length(args.reference_file)
    pseudogenome_sequence = filtered_bcf_to_fasta(args.filtered_bcf_file, reference_length)
    id = os.path.basename(args.filtered_bcf_file).split('.')[0]
    write_sequence(args.output_fasta_file, id, pseudogenome_sequence)

