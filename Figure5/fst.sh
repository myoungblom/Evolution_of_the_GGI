#!/bin/bash

# make vcf of core genome alignment using SnpSites
/opt/PepPrograms/snp-sites.2.4.1/src/snp-sites -v -o core_gene_alignment.vcf core_gene_alignment.aln

# remove ambiguous and non biallelic SNPs, reformat VCF
bash vcflibConversion.sh core_gene_alignment.vcf core

# custom script to run vcflib fst outlier analysis
python ~/scripts/scripts/scoaryToVcflibFst.py ggi_scoary_traits.csv core_gene_alignment_vcflib.vcf ggi
