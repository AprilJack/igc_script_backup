#!/bin/bash 
# run job in the current working directory where qsub is executed from
#$ -cwd
# specify that the job requires 4GB of memory
#$ -l m_mem_free=4G
#example preprocessing of seqeuncing libary ZL145

 
#unzip original datafiles                                                                                
gunzip *.gz

  
#strip fastq files and clip sequences

awk "NR%4==2" ZL198_MSeq087_SC_S1_R1_001.fastq | cut -b 1-32 > ZL198_MSeq087_SC_1_stripped.txt #barcode+YY
rm ZL198_MSeq087_SC_S1_R1_001.fastq
awk "NR%4==2" ZL198_MSeq087_SC_S1_R2_001.fastq | cut -b 1-20 > ZL198_MSeq087_SC_2_stripped.txt #12nt tag+ 8nt index
rm ZL198_MSeq087_SC_S1_R2_001.fastq



#make a new file that contains only one sequence per sequenced cluster
paste -d '' ZL198_MSeq087_SC_1_stripped.txt ZL198_MSeq087_SC_2_stripped.txt > ZL198_MSeq087_SC_PE.txt


#split dataset according to inline indexes using fastx toolkit; this by default allows up to 1 missmatch. we could go higher if we want, though maybe not neccessary
mkdir barcodesplitter
cd barcodesplitter

nl ../ZL198_MSeq087_SC_PE.txt |awk '{print ">" $1 "\n" $2}'|fastx_barcode_splitter.pl --bcfile ../barcode_v2.txt --prefix ../barcodesplitter/ --eol

#from here on, do everything for every sample individually

BCidx=($(seq 96 1 187))
for i in {1..91}; do
#filter out reads with Ns, cut off indexes and unique datafiles
awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > ZL198_MSeq087_SCprocessedBC${BCidx[$i]}.txt
#split output files into two files per index, one that is containing the read counts of each unique sequnce, the other the unique sequences themselves.
awk '{print $1}' ZL198_MSeq087_SCprocessedBC${BCidx[$i]}.txt > ZL198_MSeq087_SCBC${BCidx[$i]}_counts.txt
awk '{print $2}' ZL198_MSeq087_SCprocessedBC${BCidx[$i]}.txt > ZL198_MSeq087_SC_BC${BCidx[$i]}seq.txt
done


mkdir thresholds
cd thresholds



