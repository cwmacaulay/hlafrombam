#!/bin/bash
# HLAcounts.sh
# Charles Macaulay
# 2022-01-30

patient=$1
sampleID=$2
mapped=$3
unmapped=$4

# make a directory for this sample in /N/u/cmacaulay/20220125haplotyper/
#mkdir /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID

# symlink the files we are interested in over to this directory
#ln -s $mapped /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID
#ln -s $unmapped /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID

# grab the unmapped reads
#/N/u/cmacaulay/tools/bin/samtools view -b -f 4 $unmapped > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/unmapped.bam

# grab the reads from chromosome 6 HLA locus
#/N/u/cmacaulay/tools/bin/samtools view -b -h $mapped "chr6:28510120-33480577" > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusonly.bam

# merge these 
#/N/u/cmacaulay/tools/bin/samtools merge /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/merged.bam /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/unmapped.bam /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusonly.bam

# convert this to a fastq so we can map it to the complete HLA database as reference genome
#/N/u/cmacaulay/tools/bin/samtools fastq -0 /dev/null /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/merged.bam > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusonly.fastq

# map to the complete HLA database as reference genome

# IF IT'S RNA, need to map to the coding genome:

#if [[ "$sampleID" == *"RNA"* ]]; then
#    time /N/u/cmacaulay/tools/bwa/bwa mem /N/u/cmacaulay/20220125haplotyper/refgenome/hla_nuc.fasta /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusonly.fastq > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.sam
#else
#    time /N/u/cmacaulay/tools/bwa/bwa mem /N/u/cmacaulay/20220125haplotyper/refgenome/hla_gen.fasta /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusonly.fastq > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.sam
#fi

# convert to bam and sort
#/N/u/cmacaulay/tools/bin/samtools view -bo /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.bam /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.sam
#/N/u/cmacaulay/tools/bin/samtools sort -o /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.sorted.bam /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.bam

# generate a decoder for the HLA IDs
rm /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/*txt
/N/u/cmacaulay/tools/bin/samtools view /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/hlalocusaligned.sorted.bam | cut -f 3 | sort | uniq -c | sort -nr | awk -F " " '{print $1"\t"$2}' > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.txt
tail -n +2 /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.txt > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.txt
cat /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.txt | cut -f 2 > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAlist.txt

if [[ "$sampleID" == *"RNA"* ]]; then
    cat /N/u/cmacaulay/20220125haplotyper/refgenome/hla_nuc.fasta | grep -f /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAlist.txt | awk -F " " '{print $1"\t"$2"\t"$3}' | sed -s 's/>//g' > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAcodebreaker.txt
else
    cat /N/u/cmacaulay/20220125haplotyper/refgenome/hla_gen.fasta | grep -f /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAlist.txt | awk -F " " '{print $1"\t"$2"\t"$3}' | sed -s 's/>//g' > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAcodebreaker.txt
fi

cat /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.txt | awk '{ print $2"\t"$1 }' > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.colsreversed.txt
sort /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.colsreversed.txt > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.colsreversed.sorted.txt
sort /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAcodebreaker.txt > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAcodebreaker.sorted.txt
join /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/topHLAs.tail.colsreversed.sorted.txt /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAcodebreaker.sorted.txt | awk -F " " '{print $2"\t"$1"\t"$3"\t"$4}' | sort -nr > /N/u/cmacaulay/20220125haplotyper/$patient-$sampleID/HLAresults.txt

