#!/bin/bash 
# run job in the current working directory where qsub is executed from
#$ -cwd
# specify that the job requires 4GB of memory
#$ -l m_mem_free=4G
#example preprocessing of seqeuncing libary ZL145

cd barcodesplitter
cd thresholds
  
  
#pick thresholds from matlab plots, based on a steep drop of the 32+12 read counts. avoids too many PCR and sequnecning errors in the data themselves. note, bash seems to start counting at 0, so put a random number at the first position of this array, to get the order right.


BCidx=($(seq 96 1 187))
threshold=(0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2)

#just do a quick and dirty collpase around the UMI tags
for i in {1..91}; do
    j=${threshold[$i]} 
    echo $j
    if [ "$j" != "1" ];then
	awk '$1< '$j' {print NR}' ../ZL198_MSeq087_SCBC${BCidx[$i]}_counts.txt | head -n 1 >t
	thresh=$(cat t)
	echo $thresh
	head ../ZL198_MSeq087_SC_BC${BCidx[$i]}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ZL198_MSeq087_SC${i}quickout.txt

   else
	grep -nr ^${threshold[$i]}$  ../ZL198_MSeq087_SCBC${BCidx[$i]}_counts.txt -m 1 | cut -f1 -d":" > t
	thresh=$(cat t)
	echo $thresh

	head ../ZL198_MSeq087_SC_BC${BCidx[$i]}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ZL198_MSeq087_SC${i}quickout.txt
   fi
done


#ok, thats a lot of preprocessing done. Now we have to error correct these sequneces. The all against all mapping of barcodes in using bowtie is done in the command line. after this we move into MATLAB.


mkdir indexes

for i in {1..91}
do
echo $i
in=ZL198_MSeq087_SC${i}quickout.txt
#split off real barcodes from spike-ins

grep -v 'ATCAGTCA$' $in | grep '[TC][TC]$' > ZL198_MSeq087_SCBC${i}_quickprocessed.txt
awk '{print $1}' ZL198_MSeq087_SCBC${i}_quickprocessed.txt > ZL198_MSeq087_SC${i}_counts.txt
awk '{print $2}' ZL198_MSeq087_SCBC${i}_quickprocessed.txt > ZL198_MSeq087_SC${i}_seq.txt


nl ZL198_MSeq087_SC${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL198_MSeq087_SC_BC${i}fasta2u.txt; 
bowtie-build -q ZL198_MSeq087_SC_BC${i}fasta2u.txt indexes/BC${i}fasta2u; 
bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2u ZL198_MSeq087_SC_BC${i}fasta2u.txt bowtiealignment${i}_2u.txt
awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt;awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt
done


# now deal with spike ins

for i in {1..91}; do 
echo $i; in=ZL198_MSeq087_SC${i}quickout.txt 
grep 'ATCAGTCA$' $in > ZL198_MSeq087_SCspikes${i}_quickprocessed.txt; 
awk '{print $1}' ZL198_MSeq087_SCspikes${i}_quickprocessed.txt > ZL198_MSeq087_SCspikes${i}_counts.txt; 
awk '{print $2}' ZL198_MSeq087_SCspikes${i}_quickprocessed.txt > ZL198_MSeq087_SCspikes${i}_seq.txt;  
nl ZL198_MSeq087_SCspikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL198_MSeq087_SC_spikes${i}fasta2u.txt; 
bowtie-build -q ZL198_MSeq087_SC_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ZL198_MSeq087_SC_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt; 
awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt;
awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt; 
done




