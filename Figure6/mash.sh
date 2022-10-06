#!/bin/bash

# in the directory with all of the phage sequences (phage_seq.tar.gz) calculate mash distances
# between all sequences
for f in *.fas; do /opt/PepPrograms/Mash-2.2/mash dist $f *.fas >> mge_mash_results.txt; done
