#! /usr/bin/env python
# Given a file where the first item in the line is a sequence identifier,
# i.e. matching labels of a corresponding FASTA file, this script will
# KEEP the sequences represented in the "sequence identifiers" file and
# only write those corresponding sequences in FASTA format.


from skbio.io import read
import string
import argparse
from argparse import RawTextHelpFormatter

def parse_labels(labels_fh):
    """Read seq label IDs into dictionary.
        Incase label line has other text slice first item"""
    labels = [line.split()[0] for line in labels_fh]
    return set(labels)


def filter_seqs(fasta_ifh, fasta_ofh, seq_labels, remove_ids=False, desc=False):
    # remove seqs and keep descriptors
    if remove_ids and desc:
        for seq in fasta_ifh:
            sid = seq.metadata['id']
            if sid in seq_labels:
                continue
            else:
                seq_str = str(seq).upper()
                new_str = '>' + sid + ' ' + \
                      seq.metadata['description'] + '\n' + seq_str + '\n'
                fasta_ofh.write(new_str)
    # remove seqs and descriptors
    elif remove_ids and not desc:
        for seq in fasta_ifh:
            sid = seq.metadata['id']
            if sid in seq_labels:
                continue
            else:
                seq_str = str(seq).upper()
                new_str = '>' + sid + '\n' + seq_str + '\n'
                fasta_ofh.write(new_str)
    # keep descriptors and keep seqs
    elif desc and not remove_ids:
        for seq in fasta_ifh:
            sid = seq.metadata['id']
            if sid in seq_labels:
                seq_str = str(seq).upper()
                new_str = '>' + sid + ' ' + \
                      seq.metadata['description'] + '\n' + seq_str + '\n'
                fasta_ofh.write(new_str)
            else:
                continue
    # default: remove descriptors and keep seqs
    else:
        for seq in fasta_ifh:
            sid = seq.metadata['id']
            if sid in seq_labels:
                seq_str = str(seq).upper()
                new_str = '>' + sid + '\n' + seq_str + '\n'
                fasta_ofh.write(new_str)

def main():
    parser = argparse.ArgumentParser(
             description= 'This script will write out sequences based on \n'
             'sequence identifiers in a label file. ',
             formatter_class=RawTextHelpFormatter)
    req = parser.add_argument_group('REQUIRED')
    req.add_argument('-i', '--input_fasta', required=True, action='store',
                     help='Input fasta file.')
    req.add_argument('-l', '--input_sequence_labels', required=True, action='store',
                     help='File in which the first item in each line is'
                        ' a sequence label / identifier.')
    req.add_argument('-o', '--output_fasta', required=True, action='store',
                     help='Output fasta file.')
    optp = parser.add_argument_group('OPTIONAL')
    optp.add_argument('-d', '--include_description', action='store_true',
                      help='Boolean. Keep the additional FASTA header '
                      'description text.[Default: False]')
    optp.add_argument('-r', '--remove_ids', action='store_true',
                      help='Boolean. Remove sequences with the corresponding '
                      'IDs, rather than keep. [Default: False]')

    p = parser.parse_args()

    input_fasta = read(p.input_fasta, format='fasta')
    input_labels = open(p.input_sequence_labels, 'U')
    output_fasta = open(p.output_fasta, 'w')
    remove_ids = p.remove_ids
    include_description = p.include_description

    seq_labels = parse_labels(input_labels)
    filter_seqs(input_fasta, output_fasta, seq_labels, remove_ids=remove_ids,
               desc=include_description)

    input_fasta.close()
    output_fasta.close()

if __name__ == '__main__':
    main()
