#!/bin/bash

# For making the required input files for HMMRATAC
# usage: 
# ./Make_HMMR_Files.sh (UnsortedBAMFile) (GenomeStatsFile) (OutputPrefix)

# input: 
readFile=$1 # unsorted BAM file containing the paired end reads from an ATAC-seq experiment
genomeFile=$2 # Genome stat file with the format: ChromosomeName <TAB> ChromosomeSize
prefix=$3 #output prefix. all output files will have this prefix

# output: all in the working directory
# OutputPrefix.bam                              # Sorted BAM file of ATAC-seq reads, ready for HMMRATAC as the -b input
# OutputPrefix.bam.bai                          # Index of sorted BAM file of ATAC-seq reads, ready for HMMRATAC as the -i input
# OutputPrefix_fromBedtools.bw                  # Genome-wide bigWig file of ATAC-seq coverage, made by Bedtools, ready for HMMRATAC as the -w input
# OutputPrefix_fromMACS2.bw                     # Genome-wide bigWig file of ATAC-seq coverage, made by MACS2, ready for HMMRATAC as the -w input
# Note: the bigWig file created by MACS2 is recommended over the bedtools version. We have found that bedtools has trouble in the conversion from BAM to BED.


# executables - Please fill in with the path to each tool. 
samtools=samtools
bedGraphToBigWig=bedGraphToBigWig
bedtools=bedtools
macs=macs2

#Set parameters
sort=sort
index=index
call=callpeak
bamtobed="bamtobed -i"
genomecov=genomecov

#Sort the BAM file
#Note: this command may not work for all versions of samtools. Later versions require a '-o' for the output declaration
echo $samtools $sort $1 $3
$samtools $sort $1 $3


#Index the sorted BAM file
echo $samtools $index $3.bam $3.bam.bai
$samtools $index $3.bam $3.bam.bai


#Create Genome-wide bigWig file of ATAC-seq coverage, using bedtools, not recommended

#Use bedtools to create bigWig
echo $bedtools $bamtobed $3.bam "> TEMP_bed.bed"
$bedtools $bamtobed $3.bam > TEMP_bed.bed

echo $bedtools $genomecov "-i TEMP_bed.bed -g " $2 " > TEMP_bedgraph.bg"
$bedtools $genomecov -i TEMP_bed.bed -g  $2 -bga > TEMP_bedgraph.bg


echo $bedGraphToBigWig "TEMP_bedgraph.bg " $2 $3.bw 
$bedGraphToBigWig TEMP_bedgraph.bg $2 $3_fromBedtools.bw


echo "rm TEMP_bed.bed TEMP_bedgraph.bg"
rm TEMP_bed.bed TEMP_bedgraph.bg

#Create the Genome-wide bigWig file of ATAC-seq coverage, using MACS2, recommended


echo $macs $call -t  $3.bam --keep-dup all -n TEMP_MACS --nolambda -B
$macs $call -t $3.bam --keep-dup all -n TEMP_MACS --nolambda -B


echo $bedGraphToBigWig "TEMP_MACS_treat_pileup.bdg " $2 $3.bw
$bedGraphToBigWig TEMP_MACS_treat_pileup.bdg  $2 $3_fromMACS2.bw


echo "rm TEMP_MACS*"
rm TEMP_MACS*
