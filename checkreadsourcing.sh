#!/bin/bash
# checkreadsourcing.sh
# Charles Macaulay
# 2022-02-23

# target directory
ODIR=$1

# WHO-HLA aligned bam file
WHOHLABAM=$2

# Original hg38 bam file
ORGBAM=$3

# extract the columns from the WHO-HLA-aligned bam file for the read ID and HLAID. Reverse the columns and sort. 
/N/u/cmacaulay/tools/bin/samtools view $WHOHLABAM | cut -f 1,3 | awk '{ print $2"\t"$1}' | sort > "$ODIR"hlalocusaligned.rev.sorted.txt

# join this file with the previously-generated HLAID "codebreaker" file that translates HLAIDs into the actual 
# WHO/IMGT nomenclature so that we can pull out the reads we want. Print the output to a file that can be used
# to search the original chr6 and unmapped reads (hg38) BAM file. 
join "$ODIR"hlalocusaligned.rev.sorted.txt "$ODIR"HLAcodebreaker.sorted.txt > "$ODIR"alleles.readIDs.txt

# split the alleles.readIDs file into HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DQA1, HLA-DRA.
# HLA-1
cat "$ODIR"alleles.readIDs.txt | grep ' A\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_A.txt
cat "$ODIR"alleles.readIDs.txt | grep ' B\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_B.txt
cat "$ODIR"alleles.readIDs.txt | grep ' C\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_C.txt

# HLA-2
cat "$ODIR"alleles.readIDs.txt | grep ' DPA1\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_DPA1.txt
cat "$ODIR"alleles.readIDs.txt | grep ' DQA1\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_DQA1.txt
cat "$ODIR"alleles.readIDs.txt | grep ' DRA\*' | awk -F " " '{ print $2}' > "$ODIR"QNAMES_DRA.txt

# search the original chr6 and unmapped reads (hg38) BAM file and output results to filtered versions of the 
# original hg38 BAM files. sort & index those babies. 
# A
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_A.txt -o "$ODIR"filtered.A.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.A.bam "$ODIR"filtered.A.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.A.bam

# B
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_B.txt -o "$ODIR"filtered.B.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.B.bam "$ODIR"filtered.B.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.B.bam

# C
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_C.txt -o "$ODIR"filtered.C.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.C.bam "$ODIR"filtered.C.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.C.bam

# DPA1
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_DPA1.txt -o "$ODIR"filtered.DPA1.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.DPA1.bam "$ODIR"filtered.DPA1.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.DPA1.bam

# DQA1
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_DQA1.txt -o "$ODIR"filtered.DQA1.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.DQA1.bam "$ODIR"filtered.DQA1.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.DQA1.bam

# DRA
/N/u/cmacaulay/tools/bin/samtools view -N "$ODIR"QNAMES_DRA.txt -o "$ODIR"filtered.DRA.bam $ORGBAM
/N/u/cmacaulay/tools/bin/samtools sort -o "$ODIR"filtered.sorted.DRA.bam "$ODIR"filtered.DRA.bam
/N/u/cmacaulay/tools/bin/samtools index "$ODIR"filtered.sorted.DRA.bam

# from the filtered hg38 bam file, extract reads from A, B, C, DPA1, DQA1, DRA and count how many there are. 

# Reads from HLA-A in WHO-HLA 
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:29942532-29945870" > "$ODIR"AinA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:31353875-31357179" > "$ODIR"AinB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:31268749-31272092" > "$ODIR"AinC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:33064569-33073677" > "$ODIR"AinDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:32637406-32643684" > "$ODIR"AinDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.A.bam "chr6:32439887-32445046" > "$ODIR"AinDRA.bam

# Reads from HLA-A in WHO-HLA                                                                                                                                                                                                              
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:29942532-29945870" > "$ODIR"BinA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:31353875-31357179" > "$ODIR"BinB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:31268749-31272092" > "$ODIR"BinC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:33064569-33073677" > "$ODIR"BinDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:32637406-32643684" > "$ODIR"BinDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.B.bam "chr6:32439887-32445046" > "$ODIR"BinDRA.bam

# Reads from HLA-A in WHO-HLA                                                                                                                                                                                                              
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:29942532-29945870" > "$ODIR"CinA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:31353875-31357179" > "$ODIR"CinB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:31268749-31272092" > "$ODIR"CinC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:33064569-33073677" > "$ODIR"CinDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:32637406-32643684" > "$ODIR"CinDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.C.bam "chr6:32439887-32445046" > "$ODIR"CinDRA.bam

# Reads from HLA-A in WHO-HLA                                                                                                                                                                                                              
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:29942532-29945870" > "$ODIR"DPA1inA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:31353875-31357179" > "$ODIR"DPA1inB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:31268749-31272092" > "$ODIR"DPA1inC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:33064569-33073677" > "$ODIR"DPA1inDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:32637406-32643684" > "$ODIR"DPA1inDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DPA1.bam "chr6:32439887-32445046" > "$ODIR"DPA1inDRA.bam

# Reads from HLA-A in WHO-HLA                                                                                                                                                                                                              
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:29942532-29945870" > "$ODIR"DQA1inA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:31353875-31357179" > "$ODIR"DQA1inB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:31268749-31272092" > "$ODIR"DQA1inC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:33064569-33073677" > "$ODIR"DQA1inDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:32637406-32643684" > "$ODIR"DQA1inDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DQA1.bam "chr6:32439887-32445046" > "$ODIR"DQA1inDRA.bam

# Reads from HLA-A in WHO-HLA                                                                                                                                                                                                              
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:29942532-29945870" > "$ODIR"DRAinA.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:31353875-31357179" > "$ODIR"DRAinB.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:31268749-31272092" > "$ODIR"DRAinC.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:33064569-33073677" > "$ODIR"DRAinDPA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:32637406-32643684" > "$ODIR"DRAinDQA1.bam
/N/u/cmacaulay/tools/bin/samtools view -b -h "$ODIR"filtered.sorted.DRA.bam "chr6:32439887-32445046" > "$ODIR"DRAinDRA.bam
