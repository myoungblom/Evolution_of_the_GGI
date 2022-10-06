#!/bin/bash

# python script to remove gaps from vcf (https://github.com/tatumdmortimer/formatConverters.git)
# vcf is created using SnpSites (see fst.sh)
python gaplessVCF.py core_gene_alignment.vcf core_gene_alignment.aln

# calculating rescaling factor for TreeTime as described in the docs: https://treetime.readthedocs.io/en/latest/tutorials/homoplasy.html
# # of variable positions / length of alignment = rescaling factor
# 22839 / 1315086 = 0.01736
# set n very high so you are sure to get all homoplasies in the results
python3 /opt/PepPrograms/treetime/bin/treetime homoplasy --aln gapless.vcf --vcf-reference NCCP11945.fasta --tree GCcore_RAxML.newick --rescale 0.01736 --detailed -n 24000 > all_homoplasy_output.txt

# get the total and terminal branch legnths
grep "TOTAL tree length" all_homoplasy_output.txt
grep "TERMINAL branch length" all_homoplasy_output.txt

# mutations under the header "The ten most homoplastic mutations are:" and above the header
# "Taxons that carry positions ..." were reformatted to include
# just the mutation position and the multiplicity and written into all_gc_homoplasies.txt 
