#!/bin/bash
# try to quantify the peak density for chip with replicates, this is useful for the spike-in analysis.
# we have one flag multiple times to input more than one tag directories.
# need to redirect the final output using ">"
# 12-20-2018 update: check the repeats after getting the peaks, before annotatePeaks

ARGS=4
if [ $# -lt "$ARGS" ]
then
	echo "wrapper to quantify the peaks of ChIPs with replicates";
	echo "Usage: Quantify_ChIP_Replicates.sh -d <target tag directory1> -d <target tag directory2> -b <background tag directory1> -b <background tag directory2> -i <single input tag directory> -s <stype: histone|factor> -g <genome: hg19|mm10|rn6> > diff.out";
	exit $E_BADARGS;
fi



sample=""
control=""
input=""
style=""
genome=""


while getopts "d:b:i:s:g:" opt; do
	case $opt in
		d) sample+="$OPTARG "
		;;
		b) control+="$OPTARG "
		;;
		i) input+="$OPTARG "
		;;
		s) style="$OPTARG"
		;;
		g) genome="$OPTARG"
		;;
	esac
done 

if [ -z "$sample" ]; then
	>&2 echo "No sample info"
	exit 1
else 
	>&2 echo "sample directories:"
  >&2 echo $sample
fi
if [ -z "$control" ]; then
		>&2 echo "No control info"
		exit 1
else
	>&2 echo "control directories:"
	>&2 echo $control
fi 
if [ -z "$input" ]; then
	>&2 echo "continue with no input control"
	input_args=""
else
	>&2 echo "input directories:"
	>&2 echo $input
	input_args="-i $input"
fi
if [ -z "$style" ]; then
		>&2 echo "No style"
		exit 1
else
	>&2 echo "style:"
	>&2 echo $style
fi 
if [ -z "$genome" ]; then
		>&2 echo "No genome"
		exit 1
else
	>&2 echo "genome:"
	>&2 echo $genome
fi 

if [ "$style" = "histone" ]; then
	peakFile="regions.txt"
elif [ "$style" = "factor" ]; then
	peakFile="peaks.txt"
else
	>&2 echo "style is not histone or factor!"
	exit 1
fi 

# it should work for single input tag directory, not sure for multiple input tag directories.


# first pool the directories together 
rand_sample=$RANDOM.sample.tmp
makeTagDirectory $rand_sample -d $sample
findPeaks $rand_sample -style $style -o auto $input_args

rand_control=$RANDOM.control.tmp 
makeTagDirectory $rand_control -d $control
findPeaks $rand_control -style $style -o auto $input_args

# merge the peaks
rand_merged=$RANDOM.merged.txt
# mergePeaks -d given $rand_sample/$peakFile $rand_control/$peakFile > $rand_merged
# mergePeaks process might be a bit different from homer getDifferentialPeaksReplicates.pl, it seems that there is no -d given flag
# the two flags produced identical results
mergePeaks $rand_sample/$peakFile $rand_control/$peakFile > $rand_merged

# check repeats
repeat=/gpfs/tools/homer/data/genomes/$genome/$genome.repeats
mergePeaks -d given $rand_merged $repeat -cobound 1
# the output will be called coBoundBy0.txt and coBoundBy1.txt, and we want the coBoundBy0.txt
rand_merged_noRepeats=$RANDOM.merged.noRepeat.txt
mv coBoundBy0.txt $rand_merged_noRepeats

# get the raw count for all the samples
# annotatePeaks.pl $rand_merged  $genome -d $sample $control -noadj # it is said that -raw eq -noadj
annotatePeaks.pl $rand_merged_noRepeats  $genome -d $sample $control -raw 

rm -r $rand_sample
rm -r $rand_control
rm $rand_merged
rm $rand_merged_noRepeats


>&2 echo "done"

# Quantify_ChIP_Replicates.sh -d SB_H12ctrl_Nup93_1_S54_R1_001 -d SB_H12ctrl_Nup93_2_S60_R1_001 -d SB_H12ctrl_Nup93_3_S61_R1_001 -b SB_H12_93KD_Nup93_1_S57_R1_001 -b SB_H12_93KD_Nup93_2_S58_R1_001 -b SB_H12_93KD_Nup93_3_S59_R1_001  -s histone -g hg19 > merged_pooled_H12ctrl_and_H12_93KD_test_regions.txt