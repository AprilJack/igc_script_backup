#!/bin/sh
if [ $# -lt 1 ] 
then
	echo "Usage: script directory with fastq.gz files to use for hg19ad5 RNA-seq bedGraph generation with two different normalizations"
else
	echo "Starting hg19ad5 pipe on $1"
	runSTAR hg19ad5.star .
	cd tags
	for i in * 
	do
		makeUCSCfile $i -fragLength given -strand - -neg -o ${i}_hg19ad5_minus.bedGraph -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM chrMt
		makeUCSCfile $i -fragLength given -strand + -o ${i}_hg19ad5_plus.bedGraph -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
chrX chrY chrM chrMt
		makeTagDirectory ../ad5Only/$i -t ${i}/ad5*.tsv
		makeUCSCfile ../ad5Only/$i -fragLength given -strand - -neg -o ../ad5Only/${i}_ad5only_minus.bedGraph
		makeUCSCfile ../ad5Only/$i -fragLength given -strand + -o ../ad5Only/${i}_ad5only_plus.bedGraph
	done
	cleanSAMs .
	
fi
