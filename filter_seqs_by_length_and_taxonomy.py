#! /usr/bin/env python
# Using a minimum sequence length per taxonomic group, this script will read in
# a FASTA file and a taxonomy file, any sequence that does not fit the
# length criteria for a given taxonomic group will be discarded.
# For example, the user will specify:
#   Bacteria : 1200
#   Archaea : 900
# The means that if the sequence is from Bacteria, then it will discard any
# sequences that are less than 1200 bases long. If the sequences is from
# Archaea, then it will remove any sequence less than 900 bases long.

from skbio.io import read
from skbio import DNA
import re
import argparse
from argparse import RawTextHelpFormatter


def make_tax_group_dict(taxonomic_groups):
    """Convert string to dictionary. Then check if all values are, in fact,
    integers"""
    tg = eval(taxonomic_groups)
    if type(tg) == dict:
        pass
    else:
        raise TypeError("Taxonomic groups can not be coerced into a"
                        "dictionary! Check `--taxonomic_groups` input!")
    all_ints = all(type(value) == int for value in tg.values())
    if all_ints == True:
        return tg

def make_taxonomy_dict(taxonomy_ifh):
    d = {}
    for tax_line in taxonomy_ifh:
        id,taxonomy=tax_line.strip().split(None,1)
        d[id] = taxonomy
    return d

def filter_seqs_by_len_and_tax(fasta_ifh, fasta_ofh, taxonomic_groups_dict,
                              id_taxonomy_dict, global_length_min=1200):
    """Check if taxonomic group is present. Filter based on set sequence
    for group. Perform some taxonomy sanity checking."""
    for seq in fasta_ifh:
        seq_id_str = seq.metadata['id']
        seq_str = str(seq)

        try:
            #tax_list = id_taxonomy_dict[seq_id_str].strip().split(';')
            tax = id_taxonomy_dict[seq_id_str]
        except:
            raise KeyError("Seq ID not found in Taxonomy")
            break

        # sanity check that only one taxonomic group is present in string
        tax_list = tax.strip().split(';')
        found_group = [i for i in tax_list if i in taxonomic_groups_dict]
        lg = len(found_group)
        if lg == 0: # if no group, use global minimum seq length
            if len(seq_str) >= global_length_min:
                fasta_str = '>' + seq_id_str + '\n' + seq_str + '\n'
                fasta_ofh.write(fasta_str)
        elif lg > 1:
            print("More than one taxonomic group found in %s!" % seq_id_str )
            break
        else:
            gs = found_group[0]
            if len(seq_str) >= taxonomic_groups_dict[gs]:
                fasta_str = '>' + seq_id_str + '\n' + seq_str + '\n'
                fasta_ofh.write(fasta_str)


def main():
    parser = argparse.ArgumentParser(
             description= 'Using a minimum sequence length per taxonomic '
             'group. \nThis script will read in a FASTA file and a taxonomy file, '
             '\nany sequence that does not fit the length criteria for a given '
             '\ntaxonomic group will be discarded. For example, if the following '
             '\ncriteria are specified:\n'
             '\n\t\'{"d__Bacteria":1200, "d__Archaea":900}\'\n'
             '\nThis means, any Bacterial and Eukaryal sequences less than 1200 '
             '\nbases, and any Archaeal sequences less than 900 bases, will be '
             '\ndiscarded.',
             formatter_class=RawTextHelpFormatter)
    req = parser.add_argument_group('REQUIRED')
    req.add_argument('-i', '--input_sequences', required=True, action='store',
                     help='Input fasta file.')
    req.add_argument('-t', '--input_taxonomy', required=True, action='store',
                     help='Input taxonomy file.')
    req.add_argument('-o', '--output_sequences', required=True, action='store',
                     help='Output filtered FASTA file.')
    optp = parser.add_argument_group('OPTIONAL')
    optp.add_argument('-g', '--taxonomic_groups', action='store',
                     default='{"d__Bacteria":1200, "d__Archaea":900, "d__Eukaryota":1400}',
                     help='List of taxonomic groups and associated minimum seq '
                     '\nlength. Any sequences greater than or equal to length n.'
                     '\nTip: set to \'{}\' if you only want to use the '
                     '\n\'global_length_min\' option.'
                     "\n[Default: \'%(default)s\']")
    optp.add_argument('-m', '--global_length_min', action='store',
                     default='1200', type=int,
                     help='Any taxonomic groups not specified, will have their '
                     '\nsequences discarded if they do not fit this length '
                     '\ncritera. Set to large value if you want to remove all '
                     '\nunspecified taxonomic groups.'
                     "groups.\n[Default: %(default)s]")

    p = parser.parse_args()

    input_sequences = read(p.input_sequences, format='fasta')
    output_sequences = open(p.output_sequences, 'w')
    input_taxonomy = open(p.input_taxonomy, 'U')
    taxonomic_groups = p.taxonomic_groups
    global_length_min = p.global_length_min

    id_taxonomy_dict = make_taxonomy_dict(input_taxonomy)
    taxonomic_groups_dict = make_tax_group_dict(taxonomic_groups)

    filter_seqs_by_len_and_tax(input_sequences, output_sequences,
               taxonomic_groups_dict, id_taxonomy_dict,
               global_length_min=global_length_min)

    input_sequences.close()
    output_sequences.close()

if __name__ == '__main__':
    main()
