#!/bin/bash

#####
# This script takes in a snp vcf file produced by snp-sites
# and outputs a vcf which is properly formmated for use with vcflib,
# as well as a tabix index file. (Do not include numbers at end of referenceName!)
# Requires bgzip and tabix.
#####

# vcflibConversion.sh snpvcf.vcf referenceName

# copy original vcf to new name
p=${1/.vcf/}
vcflibFile="$p"_vcflib.vcf
cp $1 "$vcflibFile"

# separate header and contents of vcf file
head -4 "$vcflibFile" > "$p"_header.vcf
grep -vE "#" $1 > "$p"_headerless.vcf

# separate info from genotypes
cut -f1-9 "$p"_headerless.vcf > "$p"_info.vcf
cut -f10- "$p"_headerless.vcf > out.vcf

# process vcf to include diploid genotypes instead of haploid genotypes
sed 's/1\t/1\|1\t/g' out.vcf > outI.vcf
sed -z 's/\t1\n/\t1\|1\n/g' outI.vcf > outII.vcf
sed 's/0\t/0\|0\t/g' outII.vcf > outIII.vcf
sed -z 's/\t0\n/\t0\|0\n/g' outIII.vcf > outIV.vcf
sed 's/2\t/2\|2\t/g' outIV.vcf > outV.vcf
sed -z 's/\t2\n/\t2\|2\n/g' outV.vcf > outVI.vcf

# assemble new vcf
paste "$p"_info.vcf outVI.vcf > "$vcflibFile"_headerless.vcf
cat "$p"_header.vcf "$vcflibFile"_headerless.vcf > "$vcflibFile"

# remove ambiguous (*) and non-biallelic (,) sites
sed -i '/\*/d' "$vcflibFile"
#sed -i '/\,/d' "$vcflibFile"

# change reference name from "1" to referenceName
sed -i "s/^\s*1/$2/g" "$vcflibFile"

# remove temporary files
rm "$p"_header.vcf
rm "$p"_headerless.vcf
rm "$p"_info.vcf
rm out.vcf
rm outI.vcf
rm outII.vcf
rm outIII.vcf
rm outIV.vcf
rm outV.vcf
rm outVI.vcf
rm "$vcflibFile"_headerless.vcf

# make tabix index file
bgzip -c "$vcflibFile" > "$vcflibFile".gz
tabix -p vcf "$vcflibFile".gz
