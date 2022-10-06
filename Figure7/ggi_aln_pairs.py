#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter
import glob


if len(sys.argv) !=2:
	print("Usage: script.py gene-list")
	sys.exit(0)

geneList = open(sys.argv[1],"r")

###################################################
# iterate through each sequence from each gene aln
# 	for each aln in the list:
# 		pair w each aln of the list
#		write fasta (with two seq)
###################################################

def make_alns(genes):
	for index,line in enumerate(genes):
		line = line.strip().split("\t")
		gene = str(line[0])
		print(gene)
		iso_pairs = []
		alnName = "./ggi_alns/"+gene+"_aln_noStop.fasta"
		for seq1 in SeqIO.parse(alnName, "fasta"):
			iso1 = seq1.id.split("_")[0]
			for seq2 in SeqIO.parse(alnName, "fasta"):
				iso2 = seq2.id.split("_")[0]
				pairName = iso1+"-"+iso2
				revPairName = iso2+"-"+iso1
				if (iso1 != iso2) and (pairName not in iso_pairs) and (revPairName not in iso_pairs):
					tempPair = [seq1, seq2]
					tempAlnName = "./ggi_aln_pairs/"+gene+"_"+pairName+".fasta"
					with open(tempAlnName, "w"):
						SeqIO.write(tempPair, tempAlnName, "fasta")
					iso_pairs.append(pairName)
					iso_pairs.append(revPairName)
                    
make_alns(geneList)
geneList.close()
