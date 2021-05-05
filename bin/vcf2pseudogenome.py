#!/usr/bin/env python
from pysam import VariantFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
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

def calculate_reference_lengths(reference_file):
    reference_lengths = OrderedDict()
    for record in SeqIO.parse(reference_file, format = 'fasta'):
        reference_lengths[record.id] = len(record.seq)
    return reference_lengths

def calculate_gaps_to_add(gap_start_position, gap_end_position):
    return ['-'] * (gap_end_position - gap_start_position)

def filtered_bcf_to_fasta(filtered_bcf_file, reference_lengths):
    # make dictionaries to capture seeuence(s) and ongoing positons
    sequences = OrderedDict()
    for chrom in reference_lengths.keys():
        sequences[chrom] = []
    previous_positions = OrderedDict()
    for chrom in reference_lengths.keys():
        previous_positions[chrom] = 0
    with VariantFile(filtered_bcf_file) as vcf_reader:
        for record in vcf_reader.fetch():
            if record.pos % 10000 == 0:
                print(record.pos)
            record_chrom = record.chrom
            if record.pos == previous_positions[record_chrom]: # Insertion - remove previous character and add 'N'
                sequences[record_chrom].pop() # remove last position
                sequences[record_chrom].append('N')
            else:
                if previous_positions[record_chrom] + 1 < record.pos: # large deletion add gaps before adding next position
                    sequences[record_chrom].extend(calculate_gaps_to_add(previous_positions[record_chrom] + 1, record.pos))
                if 'PASS' in record.filter.keys(): # HQ SNP
                    gt = record.samples[0]['GT'][0] # get genotype (1st from tuple)
                    if 'INDEL' in record.info: # indel
                        sequences[record_chrom].append('N')
                    elif gt == 0: # reference base
                        sequences[record_chrom].append(record.ref.lower()) # add reference base as lower case
                    else: # variant
                        if len(record.alts) != 1: # if more than one ALT genotype so add N
                            sequences[record_chrom].append('N')
                        else: # add ALT SNP as upper case
                            sequences[record_chrom].append(record.alts[gt-1].upper())
                else: # if not PASS it's a low qual SNP so add N
                    sequences[record_chrom].append('N')
            previous_positions[record_chrom] = record.pos
        if previous_positions[record_chrom] != reference_lengths[record_chrom]: # if gap at the end
            sequences[record_chrom].extend(calculate_gaps_to_add(previous_positions[record_chrom], reference_lengths[record_chrom]))
        return sequences

def write_sequence(filepath, id, sequence):
    with open(filepath, 'w') as output:
        record = SeqRecord(Seq(sequence), id = id, description = '')
        SeqIO.write(record, output, "fasta")

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

    reference_length = calculate_reference_lengths(args.reference_file)
    pseudogenome_sequence_lists = filtered_bcf_to_fasta(args.filtered_bcf_file, reference_length)
    pseudogenome_sequences = [ ''.join(sequence_list) for sequence_list in pseudogenome_sequence_lists.values()]
    combined_pseudogenome_sequence = ''.join(sequence for sequence in pseudogenome_sequences)
    
    id = os.path.basename(args.filtered_bcf_file).split('.')[0]
    write_sequence(args.output_fasta_file, id, combined_pseudogenome_sequence)

