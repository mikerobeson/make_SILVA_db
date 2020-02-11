#! /usr/bin/env python
# This script will simply de-gap ('.', '-') an aligned sequence FASTA file.
# Will remove whitespace.
# Will optionally convert Us to Ts. Optionally remove the description text.
# That is, this:
# >seq1 H. Sapiens
# ..ACCGGUU---GGCCGUUCAGGGUACAGGUUGGCCGUUCAGGGUAA....
# # will be output as:
# >seq1
# ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA


from skbio.io import read
import string
import argparse
from argparse import RawTextHelpFormatter


def make_trans_table(convu=False):
    """Convert remove whitespace, gaps and optionally convert U to T."""
    replace_dict = {".":"", "-":""}
    replace_dict.update({ord(c): None for c in string.whitespace})
    if convu:
        replace_dict.update({'U':'T'})
    tt = str.maketrans(replace_dict)
    return tt

def parse_seqs(fasta_ifh, fasta_ofh, convu=False, desc=False):
    tt = make_trans_table(convu=convu)
    if desc:
        for seq in fasta_ifh:
            seq_str = str(seq).upper()
            seq_str = seq_str.translate(tt)
            new_str = '>' + seq.metadata['id'] + ' ' + \
                      seq.metadata['description'] + '\n' + seq_str + '\n'
            fasta_ofh.write(new_str)
    else:
        for seq in fasta_ifh:
            seq_str = str(seq).upper()
            seq_str = seq_str.translate(tt)
            new_str = '>' + seq.metadata['id'] + '\n' + seq_str + '\n'
            fasta_ofh.write(new_str)


def main():
    parser = argparse.ArgumentParser(
             description= 'This script will simply degap FASTA files.\n'
             'Optionally without the description and or converting '
             'Us to Ts.\n'
             'That is, this: \n'
             '\t>seq1 H. Sapiens\n'
             '\t...ACCGGUU---GGCCGUU CAGGGUACAGGUUGGCCGUUCAGGGUAA...\n'
             'will be output as:\n'
             '\t>seq1\n'
             '\tACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA\n',
             formatter_class=RawTextHelpFormatter)
    req = parser.add_argument_group('REQUIRED')
    req.add_argument('-i', '--input_fasta', required=True, action='store',
                     help='Input fasta file.')
    req.add_argument('-o', '--output_fasta', required=True, action='store',
                     help='Output fasta file.')
    optp = parser.add_argument_group('OPTIONAL')
    optp.add_argument('-d', '--include_description', action='store_true',
                      help='Boolean. Keep the additional FASTA header '
                      'description text.[Default: False]')
    optp.add_argument('-u', '--convert_to_uracil', action='store_true',
                      help='Boolean. Convert "U" to "T". [Default: False]')

    p = parser.parse_args()

    input_fasta = read(p.input_fasta, format='fasta')
    output_fasta = open(p.output_fasta, 'w')
    convert_to_uracil = p.convert_to_uracil
    include_description = p.include_description

    parse_seqs(input_fasta, output_fasta, convu=convert_to_uracil,
               desc=include_description)

    input_fasta.close()
    output_fasta.close()

if __name__ == '__main__':
    main()
