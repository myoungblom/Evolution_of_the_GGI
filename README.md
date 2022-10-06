# Evolution of the GGI
Scripts and data associated with Figures 5-8 of GGI manuscript.
Working title: "The Gonococcal Genetic Island defines distinct sub-populations of *Neisseria gonorrhoeae*"  

### Figure 5
**core\_gene\_alignment[...]**: *N. gonorrhoeae* core genome alignment and related files including VCF from SnpSites and VCF converted for use with Vcflib by vcflibConversion.sh  
**ggi\_scoary\_traits.csv**: ggi presence absence  
**scoaryToVcflibFst.py**: runs vcflib outlier analysis and calculates null distribution using phenotype file and vcf   
**vcflibConversion.sh**: converts SnpSites VCF for use with Vcflib  
**fst.sh** - steps for running fst outlier analysis  
**GCFstManhattan.R**: script for making Figure 5 using fst oulier and homoplasy results  

##### homoplasy/
**GCCore\_RAxML.newick**: core genome phylogeny produced using RAxML  
**gappless.vcf**: vcf with gaps removed using script from https://github.com/tatumdmortimer/formatConverters.git  
**NCCP11945.fasta**: reference strain genome  
**treetime.sh**: running treetime to calculate homoplasies and tree branch lengths  

### Figure 6
**phage_seq.tar.gz**: all phage sequences as annotated by ProphET  
**mash.sh**: calculating mash distances  
**mash\_results.txt**: mash distances for all phage sequences  
**mge\_annot.tar.gz**: MGE as annotated by MobileElementFinder  
**mge\_summary.txt**: summary of number of IS elements per strain  
**GC\_phageMDS\_MGE.R**: R script for plotting Figure 6 using phage MDS and MGE results  

### Figure 7
**piNpiS.sh**: bash script for running piN/piS analyses  
**ggi\_alns.tar.gz**: GGI gene alignments  
**ggi\_gene\_order.txt**: order of GGI genes based on strain MS11  
**ggi\_aln\_pairs.py**: python script which makes fasta files for each pairwise comparison within each GGI gene  
**GGI\_piNpiS.txt**: piN/piS results for GGI genes  
**selectionStats\_piNpiS.py**: python script that calculates piN/piS using Egglib  
**GGI\_piNpiS.R**: R script for plotting piN/piS values (Figure 7C)  

### Figure 8
**CoreGGI\_aln.fasta**: core GGI genome alignment  
**GCcore\_RAxML.newick**: core genome (for *N. gonorrhoeae* not the GGI core) phylogeny produced using RAxML  
**GGIcore\_RAxML.newick**: GGI core genome phylogeny produced using RAxML  
**CoreGGI\_AncRecon.R**: R script for producing analyses and plots for Figure 8  
