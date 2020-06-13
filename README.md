# Note:
Much of this code has been refactored, and much improved, as part of the [RESCRIPt plugin](https://github.com/bokulich-lab/RESCRIPt) for QIIME 2 (Released on June 12, 2020). Give it a try! I will leave this repo here for a while. 


# make_SILVA_db
General procedure for making QIIME 2 compatible SILVA reference files

# SILVA-dbs

This repository is intended to be a collection of formatted [SILVA](https://www.arb-silva.de/) files for use in [QIIME 1](http://qiime.org/) or [QIIME 2](https://qiime2.org/). If you use the SILVA reference files be sure to read their [license](https://www.arb-silva.de/silva-license-information). The approach I take here is partly inspired by my prior experiences parsing reference databases as well as from many other discussions and online resources.


In brief, I use the following collection of scripts, with modification, and some manual edits to generate a QIIME formatted SSU & LSU gene sequence reference sets. If you are here, I assume you have QIIME 1 & QIIME 2 installed. We'll be making use of these two environments as they contain most the dependencies you need to run the pipeline outlined here.

I will eventually post more details of the pipeline and code used to generate these SSU and LSU files as time permits. Please consider the below my scratch notes.

## SILVA SSU Processing
*These steps should also work with LSU files too.*


1. Downloaded the following files from the [SILVA v138 release](https://www.arb-silva.de/no_cache/download/archive/release_138 ). Note, To save time, and avoid reclustering myself, I decided to leverage the SILVA SSU NR99 reference as outlined [here](https://www.arb-silva.de/projects/ssu-ref-nr/)

  * `tax_slv_ssu_138.txt`
  * `tax_slv_ssu_138.tre`
  * `taxmap_slv_ssu_ref_nr_138.txt`
  * `SILVA_138_SSURef_NR99_tax_silva_trunc.fasta`
  * `SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta`


2. Parse the SILVA 138 Taxonomy. (In QIIME 2 environment.)

  ```
  python parse_silva_taxonomy.py \
    -t tax_slv_ssu_138.txt \
    -p tax_slv_ssu_138.tre \
    -m taxmap_slv_ssu_ref_nr_138.txt \
    -s \
    -o SILVA_138_Taxonomy.txt
```

  * Note, I've used the optional flag to include the species labels. But be wary! For example, there are taxa annotated with the species label of the host or source rather than the sequence it self. Here is an example:
    - d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Serratia; s__Oryza_sativa
    - d__Eukaryota; p__Arthropoda; c__Insecta; o__Hemiptera; f__Hemiptera; g__Hemiptera; s__Oryza_sativa

  As you can see, we have an insect and a bacterial sequence (*Note: for the bacteria this is not a mitochondria / chloroplast / plastid sequence!*) both annotated with the species label *Oryza sativa* (rice). In most cases the species rank information seems okay, but there are enough issues like the one above, that convinced me to generally be cautious of the species label. If you do not want to make use of the species labels simply remove the `-s` flag.

  Also note, I've not curated the species names. This is important as you may have (nearly) identical sequences that point to very slightly different species label annotations, such as:
  - s__Clostridioides_difficile
  - s__Clostridioides_difficile_R20291

  So, if your sequence is similar to these, you'd think it should be classified as `s__Clostridioides_difficile`. This will not be the case, as the specific species strings are different. What the classifier may actually return is the upper-level taxonomy `g__Clostridioides`. This is not the fault of the classifier *per se*, but a problem of annotation. Because of this I only return the first two words (i.e. *Clostridioides* and *difficile*) of the "species" string.


3. Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA. Do this for both the aligned and unaligned FASTA files. (In QIIME 2 environment.)

  ```
  python convert_rna_to_dna.py \
    -i SILVA_138_SSURef_NR99_tax_silva_trunc.fasta \
    -o SILVA_seqs.fasta
  ```

  For alignment files, use `-g` to convert '.' to '-'
  ```
  python convert_rna_to_dna.py \
    -i SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta \
    -o SILVA_align_seqs.fasta \
    -g
  ```


4. Remove sequences that contain 5 or more ambiguous bases and/or homopolymers with 8 or more bases. See help text for more info. (In QIIME 2 environment.)

  ```
  python remove_seqs_with_homopolymers.py \
    -i SILVA_seqs.fasta \
    -o SILVA_seqs_polyfilt.fasta
  ```


5. Filter sequences by length, based on taxonomy. Rather than blindly filter all of the reference sequences below a certain length we'll differentially filter based on the taxonomy of the reference sequence. That is, if we decide to remove any sequences below 1200 or 1400 bp, then many of the Archaea will be lost. However, we also do not want to increase the retention of bad Bacterial or Eukaryal sequences. So, we'll mitigate these issues by differentially filtering based on length. By default the script will remove Archaeal (16S) sequences less than 900 bp, Bacterial (16S) sequences less than 1200 bp, and any Eukaryal (18S) sequences less than 1400 bp. See help text for more info. (In QIIME 2 environment.)


```
  python filter_seqs_by_length_and_taxonomy.py \
    -i SILVA_seqs_polyfilt.fasta \
    -o SILVA_seqs_polyfilt_lenfilt.fasta \
    -t SILVA_138_Taxonomy.txt
```


6. Now we'll parse the aligned FASTA file so that it contains the same sequences of our unaligned filtered sequences.


  ```
  egrep '^>' SILVA_seqs_polyfilt_lenfilt.fasta | sed 's/>//g' > ids_to_keep.txt

  python filter_fasta_by_seq_id.py \
    -l ids_to_keep.txt \
    -i SILVA_align_seqs.fasta \
    -o SILVA_align_seqs_polyfilt_lenfilt.fasta
  ```


7. Extract V4 region using [EMP 515-806 primers](http://www.earthmicrobiome.org/protocols-and-standards/16s/) locations. This approach allows us to retain more sequences within this region as opposed to using primer sequence to find and remove the corresponding region. (In QIIME 2 environment.)

  * We'll do this by making a temporary small alignment file to map primers to:

  ```
  head -n 500 SILVA_align_seqs_polyfilt_lenfilt.fasta > short_alignment.fasta
  ```

  * Run mafft to map our primers to the short alignment.

  ```
  mafft \
  --addfragments emp_primers.fasta  \
  --mapout short_alignment.fasta \
  > emp_primers_aln_map.txt
  ```

  * Let's look at the output. *Note the mapping file name is based on the `emp_primers.fasta` file name, not the output of the actual alignment*.

  ```
  cat /Users/robesonmichael/Documents/tmp/SILVA_138_Files/emp_primers.fasta.map
	>Forward_515_GTGYCAGCMGCCGCGGTAA
	# letter, position in the original sequence, position in the reference alignment
	g, 1, 11895
	t, 2, 11897
	g, 3, 13125
	y, 4, 13127
	c, 5, 13129
	a, 6, 13130
	g, 7, 13133
	c, 8, 13135
	m, 9, 13139
	g, 10, 13142
	c, 11, 13144
	c, 12, 13147
	g, 13, 13148
	c, 14, 13152
	g, 15, 13854
	g, 16, 13855
	t, 17, 13858
	a, 18, 13859
	a, 19, 13861
	>Reverse_806_revcomp_GGACTACNVGGGTWTCTAAT
	# letter, position in the original sequence, position in the reference alignment
	a, 1, 23446
	t, 2, 23447
	t, 3, 23958
	a, 4, 23959
	g, 5, 23961
	a, 6, 23963
	w, 7, 23964
	a, 8, 23965
	c, 9, 25277
	c, 10, 25283
	c, 11, 25285
	b, 12, 25287
	n, 13, 25290
	g, 14, 25292
	t, 15, 25293
	a, 16, 25294
	g, 17, 25298
	t, 18, 25300
	c, 19, 25316
	c, 20, 25318
  ```

  You can sanity-chack the alignment in your favorite alignment viewer.

  * Anyway, this looks good, lets extract amplicon region from the original alignment. *Note the values I use for the start and end position of the alignment*.

  ```
  python extract_alignment_region.py \
    -i SILVA_align_seqs_polyfilt_lenfilt.fasta \
    -o SILVA_align_seqs_polyfilt_lenfilt_empv4.fasta \
    -s 13862 \
    -e 23445
  ```

  * Remove gaps from sequence. We will be using this output to make a EMP V4 region classifier.

  ```
  python degap_fasta.py \
    -i SILVA_align_seqs_polyfilt_lenfilt_empv4.fasta \
    -o SILVA_empv4.fasta
  ```

  * Because we use the alignment to extract the V4 region, not all sequences will have enough data through this region of the alignment and must be removed; i.e. they contian only gaps). We'll keep any sequence with at least 200 bp using
  [vsearch](https://github.com/torognes/vsearch).

  ```
vsearch --fastx_filter SILVA_empv4.fasta \
  --fastq_minlen 200 \
  --fastaout SILVA_empv4_emptyrem.fasta
  ```

  * OPTIONAL STEP: Depending on how you wish to leverage your taxonomy, you may or may not want to dereplicate your sequences. If you do dereplicate your sequences, then you'll have to make a consensus taxonomy. This is a good move if you are using closed-reference style approaches, as you only want a single taxonomy for a given unique sequence variant. Otherwise you can keep the identical sequences and let the classifier, and tools like [clawback](https://github.com/BenKaehler/q2-clawback/blob/master/README.md) handle the case for identical short sequences with differing taxonomies.

    Dereplicate amplicon region

    ```
    vsearch \
      --derep_fulllength SILVA_empv4_emptyrem.fasta \
      --output SILVA_empv4_emptyrem_derep.fasta \
      --uc SILVA_empv4_emptyrem_derep.uc \
      --threads 4 \
      --fasta_width 0
    ```

    Switch to QIIME 1 environment. Use the [parse_otu_mapping_from_uc.py](https://gist.github.com/walterst/8b88b149a08ef91651f85b088efda1e2) and [create_consensus_taxonomy.py](https://gist.github.com/walterst/bd69a19e75748f79efeb) code from Tony Walters.

    ```
    python parse_otu_mapping_from_uc.py \
          SILVA_empv4_emptyrem_derep.uc \
          SILVA_empv4_emptyrem_derep_otu_map.txt


    python create_consensus_taxonomy.py \
      SILVA_138_Taxonomy.txt  \
      SILVA_empv4_emptyrem_derep.fasta  \
      SILVA_empv4_emptyrem_derep_otu_map.txt \
      SILVA_empv4_consensus_taxonomy.txt
      ```

  You can also try [create_majority_taxonomy.py](https://gist.github.com/walterst/f6f08f6583bb320bb10d).


8. Switch back to QIIME 2 environment and import files.
  * Import consensus taxonomy file for V4 region and full-length sequences:

  ```
  qiime tools import \
    --type FeatureData[Taxonomy] \
    --input-path SILVA_empv4_consensus_taxonomy.txt \
    --input-format HeaderlessTSVTaxonomyFormat \
    --output-path SILVA-v138-515f-806r-taxonomy.qza

   qiime tools import \
    --type FeatureData[Taxonomy] \
    --input-path SILVA_138_Taxonomy.txt  \
    --input-format HeaderlessTSVTaxonomyFormat \
    --output-path Silva-v138-full-length-seq-taxonomy.qza
```
  * Import FASTA files:

  ```
  qiime tools import \
    --input-path SILVA_seqs_polyfilt_lenfilt.fasta \
    --output-path SILVA-138-SSURef-Full-Seqs.qza \
    --type 'FeatureData[Sequence]'

  qiime tools import \
    --input-path SILVA_empv4_emptyrem_derep.fasta \
    --output-path SILVA-138-SSURef-515f-806r-Seqs.qza \
    --type 'FeatureData[Sequence]'
  ```


9. Train classifiers for V4 and full-length.

  ```
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads SILVA-138-SSURef-Full-Seqs.qza \
    --i-reference-taxonomy Silva-v138-full-length-seq-taxonomy.qza \
    --o-classifier SILVA-138-SSURef-full-length-classifier.qza

  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads SILVA-138-SSURef-515f-806r-Seqs.qza \
    --i-reference-taxonomy SILVA-v138-515f-806r-taxonomy.qza \
    --o-classifier SILVA-v138-515f-806r-classifier.qza
```

That's it! I'll periodically update this as time permits. If you have suggestions for improving this, then by all means, let me know. Even better, submit a pull request! If you find this useful, I'd appreciate any acknowledgments.
