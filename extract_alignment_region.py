#! /usr/bin/env python
# This script will extract a region of an aligment

from skbio.io import read as read_fasta
import argparse



def extract_region(seq_str, startp, endp):
	return seq_str[int(startp):int(endp)]


def iter_seqs(inf, outf, startp, endp):
	fh = read_fasta(inf, format='fasta')
	ofh = open(outf, 'w')

	for seq in fh:
		header = seq.metadata['id'] + ' ' +seq.metadata['description']
		extract_seq = extract_region(str(seq), startp, endp)
		new_fasta_str = '>' + header + '\n' + extract_seq + '\n'
		ofh.write(new_fasta_str)

	fh.close()
	ofh.close()


def main():
	parser = argparse.ArgumentParser(
			 description= 'Extracts region from an alignment by'
			  			   'column position.')

	req = parser.add_argument_group('REQUIRED')
	req.add_argument('-i', '--input_alignment', required=True, action='store',
	 			     help='Input alignment file.')
	req.add_argument('-o', '--output_alignment', required=True, action='store',
				     help='Output extracted alignment file')
	req.add_argument('-s', '--start_position', required=True, action='store',
	 			     type=int, help='Starting alignment column position')
	req.add_argument('-e', '--end_position', required=True, action='store',
					 type=int, help='Ending alignment column position')


	p = parser.parse_args()

	input_alignment = p.input_alignment
	output_alignment = p.output_alignment
	start_position = p.start_position-1
	end_position = p.end_position

	iter_seqs(input_alignment, output_alignment, start_position,
			  end_position)


if __name__ == '__main__':
        main()
