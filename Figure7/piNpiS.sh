#!/bin/bash

# make output directory for alignment pairs
mkdir ggi_aln_pairs

# run python script to create pairwise alignments for all ggi genes
python ggi_aln_pairs.py ggi_gene_order.txt

# run python script which uses egglib to calculate piN/piS for all pairwise alignments
python3 selectionStats_piNpiS.py -d ggi_aln_pairs
