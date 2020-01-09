#!/ur/bin/env python
# By: Mike Robeson Dec 20, 2019
# I ran this code within the `qiime2-2019.10` environment.
# Simple concept code to prepare a Greengenes-like taxonomy for SILVA (v138).

from skbio.tree import TreeNode
import re
import argparse

#allowed_ranks_list = [('domain','d__'), ('kingdom','k__'), ('phylum','p__'),
#					  ('class','c__'), ('order','o__'), ('family','f__'),
#					  ('genus','g__')]
allowed_ranks_list = [('domain','d__'), ('phylum','p__'),
					  ('class','c__'), ('order','o__'), ('family','f__'),
					  ('genus','g__')]
allowed_ranks_dict = dict(allowed_ranks_list)
allowed_ranks = allowed_ranks_dict.keys()
ranks = [ranktax[0] for ranktax in allowed_ranks_list]
rank_prefixes = [ranktax[1] for ranktax in allowed_ranks_list]

whitespace_pattern = re.compile(r'\s+')
allowed_chars = set('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-[]()/.\\')
#odd_chars = set(["'","{","}","[","]","(",")","_","-","+","=","*","&","^","%",
#				 "$","#","@","\"","/","|","`","~",':',';',",",".","?"])

# make taxonomy ID dictionary
def make_taxid_dict(taxonomy_file):
	"""Returns the dict: {TaxonomyID : (TaxonomyRank, Taxonomy)}
						  {"3698" : ("genus", "Haemophilus") }
	"""
	d = {}
	for line in taxonomy_file:
		sline = line.strip()
		if sline == '':
			continue
		else:
			tax,tid,rank = sline.split('\t')[0:3]
			tax = tax.strip()
			tid = tid.strip()
			rank = rank.strip()
			rtax = tax.rsplit(';')[-2]
			d[tid] = (rank,rtax)
	print('Number of taxonomy IDs: ', len(d))
	return d

# make accession to taxonomy id dict
def make_acc_to_species_tid_dict(taxmap_file):
	"""Returns the dict: {FullAccession : (Species, TaxonomyID)}
						{"A16379.1.1485" : ("[Haemophilus] ducreyi", "3698")}
	"""
	acc_species_tid = {}

	for line in taxmap_file:

		sline = line.strip()

		if sline.startswith("primaryAccession"):
			continue
		elif sline == '':
			continue
		else:
			ll = sline.split('\t')
			full_acc = '.'.join(ll[0:3])
			species = ll[4].strip()
			tid = ll[5].strip()
			acc_species_tid[full_acc] = (species, tid)
	print('Number of \"{Full Accession: (species, TaxID)}\" records in taxmp: ',
	      len(acc_species_tid))
	return acc_species_tid


def filter_characters(lin_name, allowed_chars=allowed_chars,
					 whitespace_pattern=whitespace_pattern):
	""" Only keep allowed characters. Should remove funny ascii too.
	Partial idea taken from https://gist.github.com/walterst/0a4d36dbb20c54eeb952
	WARNING: may result in lineage names missing characters"""

	updated_lineage_name = ""
	for char in lin_name.strip():
		if char in allowed_chars or char.isspace():
			updated_lineage_name += char

	new_lin_name = whitespace_pattern.sub("_", updated_lineage_name.strip())

	return new_lin_name


def build_base_silva_taxonomy(tree_file, tax_dict):
	"""Returns {TaxonomyID : [(rank, taxonomy), ...]} """
	print("Building base SILVA taxonomy...")
	tree = TreeNode.read(tree_file)
	ml = {}
	for node in tree.postorder():# tree.tips():
		if node.is_root():
			break

		l = []
		rank, taxonomy = tax_dict[node.name]
		clean_taxonomy_str = filter_characters(taxonomy)

		if rank in allowed_ranks:
			l.append((allowed_ranks_dict[rank], clean_taxonomy_str))

		for ancestor in node.ancestors():
			if ancestor.is_root():
				break
			else:
				arank, ataxonomy = tax_dict[ancestor.name]
				cleaned_ataxonomy = filter_characters(ataxonomy)
				if arank in allowed_ranks:
					l.append((allowed_ranks_dict[arank], cleaned_ataxonomy))

		#l.reverse()
		ml[node.name.strip()] = dict(l)

	return ml

def propagate_upper_taxonomy(sts_dict, rank_prefixes):
	print('Propagating upper level taxonomy...')
	prop_tax_dict = {}
	curr_tax = 'NOAVAILABLETAXONOMY'
	for tid, rank_taxonomy in sts_dict.items():
		prop_ranks = [''] * len(rank_prefixes)
		for i,rp in enumerate(rank_prefixes):
			try:
				tax = rank_taxonomy[rp]
				curr_tax = tax
			except:
				tax = curr_tax
			prop_ranks[i] = rp + tax
		prop_tax_dict[tid] = '; '.join(prop_ranks)
	return prop_tax_dict

def write_tax_strings(facc_species_tid_dict, prop_dict, outfile,
					  sp_label=False):
		print('Saving new fixed-rank SILVA taxonomy to file...')
		if sp_label:
			for facc,taxinfo in facc_species_tid_dict.items():
				tp = prop_dict[taxinfo[1]]
				species_name = taxinfo[0]
				clean_species_name = '; s__' + filter_characters(species_name)
				nts = facc + '\t' + tp + clean_species_name + '\n'
				#print(nts)
				outfile.write(nts)
		else:
			for facc,taxinfo in facc_species_tid_dict.items():
				tp = prop_dict[taxinfo[1]]
				nts = facc + '\t' + tp  + '\n'
				#print(nts)
				outfile.write(nts)


def main():
	parser = argparse.ArgumentParser(
	         description='Creates a GreenGenes-like formatted taxonomy for SILVA.')
	req = parser.add_argument_group('REQUIRED')
	req.add_argument('-t', '--taxonomy', required=True, action='store',
	                help='SILVA taxonomy to Accession ID file')
	req.add_argument('-p', '--taxonomy_tree', required=True, action='store',
	                 help='SILVA taxonomic hierarchy file')
	req.add_argument('-m', '--taxonomy_map', required=True, action='store',
	                 help='SILVA Accession taxonomy map')
	req.add_argument('-o', '--output_taxonomy', required=True, action='store',
	                 help='Newly formatted SILVA taxonomy')
	opt = parser.add_argument_group('OPTIONAL')
	opt.add_argument('-s', '--include_species', action='store_true',
	                    help='Boolean. Use to extract and append Species '
	                         'labels to the formatted taxonomy. WARNING: '
	                         'Species labels may not be accurate! '
							 '[Default: False]')

	#parser.print_help()
	#parser.parse_args([])

	p = parser.parse_args()

	input_taxonomy = open(p.taxonomy, 'U')
	input_taxonomy_tree = open(p.taxonomy_tree, 'U')
	input_taxonomy_map = open(p.taxonomy_map, 'U')
	sp_label = p.include_species
	ouput_taxonomy = open(p.output_taxonomy, 'w')


	tax_dict = make_taxid_dict(input_taxonomy)
	sts_dict = build_base_silva_taxonomy(input_taxonomy_tree,
										 tax_dict)
	prop_dict = propagate_upper_taxonomy(sts_dict, rank_prefixes)
	taxmap_dict = make_acc_to_species_tid_dict(input_taxonomy_map)

	write_tax_strings(taxmap_dict, prop_dict, ouput_taxonomy, sp_label=sp_label)
	ouput_taxonomy.close()



if __name__ == '__main__':
	main()
	# Example, obtain the files below from:
	# https://www.arb-silva.de/no_cache/download/archive/release_132/Exports/taxonomy/
	#	--taxonomy = "tax_slv_ssu_132.txt"
	#	--taxonomy_map = "taxmap_slv_ssu_ref_132-corrected.txt"
	#	--taxonomy_tree = "tax_slv_ssu_132.tre"
	#	There is an accession to taxonomy id file = "tax_slv_ssu_132.acc_taxid"
	#		I did not make use of this file as the Accessions and Taxonomy
	#		IDs are in the taxmap_slv* file. Which we have to parse anyway, if
	#       we want the Species labels.
