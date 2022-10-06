#!/usr/bin/env python

import sys
import os
import itertools
import subprocess
import numpy as np

#####
# This script takes in a scoary traits file and uses it to 
# bin isolates into 'target' and 'background' for use with vcflib-wcFst.
# Will also permute data to create null distribution from which
# Fst outliers can be identified.
# Requirements: vcf file must be formmated for vcflib,
# must have associated tabix index file (see vcflibConversion.sh) 
# NOTE: path to vcflib executable has been hardcoded (lines 70 & 88), change to the path in your filesystem.
#####

# check for correct commandline arguments
if len(sys.argv) != 4:
    print("Usage: scoaryToVcflibFst.py <scoaryTraits.csv> <snpvcf.vcf> <headerToParse>")
    sys.exit(0)

scoary = sys.argv[1]
vcfFile = sys.argv[2]
header = str(sys.argv[3])

backgroundList = []
targetList = []
phenotypes = {}
indexes = {}

# make dictionary with binary phenotypes from scoary traits file
total = 0
with open(scoary, 'r') as scoaryFile:
    firstline = next(scoaryFile)
    indexed = firstline.strip().split(",").index(header)
    for line in scoaryFile:
        parser = line.strip("\n").strip("\r").split(",")
        total += 1
        value = parser[indexed]
        isolate = parser[0]
        phenotypes[isolate] = value

# make dictionary with indices from vcf file
with open(vcfFile, 'r') as vcf:
    for line in vcf:
        line = line.strip("\n")
        if line.startswith("#CHROM"):
            parser = line.split("\t")
            for item in parser[9:]:
                indexes[item] = (parser.index(item))

# read from phenotypes dictionary, write from indexes dictionary into background
# and target lists
for key, value in phenotypes.items():
    if value == "0":
        backgroundList.append(str(indexes[key]-9))
    elif value == "1":
        targetList.append(str(indexes[key]-9))
targetNum = len(targetList)

# join lists to make csv
backgroundArg = ",".join(backgroundList)
targetArg = ",".join(targetList)

# use background and target lists to run vcflib wcFst
outfile = header+"_wcFst.txt"
with open(outfile, 'w') as fstOutput:
    subprocess.call(["/opt/PepPrograms/vcflib/bin/wcFst" ,"--target",targetArg,"--background",\
        backgroundArg,"--file",vcfFile,"--type", "GT"], stdout=fstOutput)

# make list of vcf indices
isolateList = []
for i in range(total):
    isolateList.append(i)

# run 100 Fst permutations based on target size
for i in range(100):
    nullOutput = open(str(i)+"_wcFst.txt", 'w')
    target = list(np.random.choice(isolateList, size=targetNum, replace=False))
    background = []
    for iso in isolateList:
        if iso not in target:
           background.append(iso)
    targetNull = ",".join(str(x) for x in list(target))
    backgroundNull = ",".join(str(x) for x in list(background))
    subprocess.call(["/opt/PepPrograms/vcflib/bin/wcFst","--target",targetNull,"--background",\
            backgroundNull, "--file", vcfFile, "--type", "GT"], stdout=nullOutput)
    nullOutput.close()

# find maximum wcFst value in null distribution
allFst = []
for i in range(100):
    with open(str(i)+"_wcFst.txt", "r") as f:
        for line in f:
            line = line.strip("\n")
            info = line.split("\t")
            allFst.append(float(str(info[4])))

# write summary file for null distribution
allFst.sort()
nullMin = min(allFst)
nullMax = max(allFst)

with open(header+"_nullWcFst.txt",'w') as nullOut:
    nullOut.write("Min\tMax\n")
    nullOut.write("%s\t%s\n" % (nullMin, nullMax))

# remove intermediate files
for i in range(100):
    if os.path.exists(str(i)+"_wcFst.txt"):
        os.remove(str(i)+"_wcFst.txt")
