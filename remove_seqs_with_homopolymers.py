#! /usr/bin/env python
# This script will read in a FASTA file and remove any sequences
# that do not meet the filtering criteria:
#   - Exsessive homopolymers
#   - Excessive ambiguous bases


from skbio.io import read
from skbio import DNA
import re
import argparse
from argparse import RawTextHelpFormatter

def filter_seqs_with_ambiguous_bases(seq, n_ambiguous_bases):
    s = DNA(seq)
    ambig_bases_in_seq = sum(s.degenerates())
    if ambig_bases_in_seq >= n_ambiguous_bases:
        return True
    else:
        return False

def filter_homopolymer(seq, n_homopolymer_length):
    nhl = n_homopolymer_length - 1 # due to how regex is written
    if nhl < 1:
        raise ValueError("Homopolymer length must be >= 2!")
    else:
        regex_str = "([ACGTURYSWKMBDHVN])\\1{%s,}" % nhl
        for p in re.finditer(regex_str, seq):
            if len(p.group()) >= n_homopolymer_length:
                return True
            else:
                continue
        return False

def filter_seqs(fasta_ifh, fasta_ofh, n_homopolymer_length=8,
                n_ambiguous_bases=5):
    for seq in fasta_ifh:
        seq_str = str(seq)
        ambig = filter_seqs_with_ambiguous_bases(seq_str, n_ambiguous_bases)
        if ambig == False:
            poly = filter_homopolymer(seq_str, n_homopolymer_length)
            if poly == False: # if we make it here, write seq to file
                seq_str = '>' + seq.metadata['id'] + ' ' + seq.metadata['description'] + '\n' + seq_str + '\n'
                fasta_ofh.write(seq_str)
        #    else:
        #        continue
        #else:
        #    continue

def main():
    parser = argparse.ArgumentParser(
             description= 'This script will read in a FASTA file and remove '
              'any sequences that have homopolymers and ambiguous base calls.'
              'That is, the following sequences would be removed: \n'
              '\t>seq1-homopolymeric\n'
              '\tACCGGTTGGCCGTTTTTTTTTCAGGGMACAGGTTVGCCGTTCAGGGTAA\n'
              '\t>seq2-ambiguos-bases\n'
              '\tACCGGTTGGCCVTGCCGMMTTCVVAGRGTAY\n',
              formatter_class=RawTextHelpFormatter)
    req = parser.add_argument_group('REQUIRED')
    req.add_argument('-i', '--input_fasta', required=True, action='store',
                     help='Input fasta file.')
    req.add_argument('-o', '--output_fasta', required=True, action='store',
                     help='Output fasta file.')
    optp = parser.add_argument_group('OPTIONAL')
    optp.add_argument('-p', '--n_homopolymer_length', action='store',
                     type=int, default=8,
                     help='Remove sequences that contain homopolymers of '
                     'greater than or equal to length n. \n'
                     "[Default %(default)s)]")
    optp.add_argument('-a', '--n_ambiguous_bases', action='store', type=int,
                     default=5, help='Remove sequences that contain a '
                     'number of IUPAC ambiguous bases greater than or equal '
                     "to length n. \n[Default %(default)s)]")

    p = parser.parse_args()

    input_fasta = read(p.input_fasta, format='fasta')
    output_fasta = open(p.output_fasta, 'w')
    n_homopolymer_length = p.n_homopolymer_length
    n_ambiguous_bases = p.n_ambiguous_bases

    filter_seqs(input_fasta, output_fasta,
                n_homopolymer_length=n_homopolymer_length,
                n_ambiguous_bases=n_ambiguous_bases)

    input_fasta.close()
    output_fasta.close()

if __name__ == '__main__':
    main()
