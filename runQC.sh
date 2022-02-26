#!/bin/bash
if [ $# -lt 1 ] 
then
	echo "runQC dirWithFastqs1 [dirWithFastqs2] ..." 
	echo "Will randomly subsample 100k reads and run fastqc, fastqScreen, and fastqBlast and output to QC"
else
for d in $@; do
  fastqfiles=$(find ${d} -name *fastq.gz)
  for f in $fastqfiles; do #for each fastq file
    dir=${f%/*}
    echo "Dir is $dir"
    rm -r $dir/QC
    mkdir $dir/QC
    sampled=${f%.fastq.gz}_sampled.fastq
    data=${sampled%.fastq}_fastqc/fastqc_data.txt
    zip=${sampled%.fastq}_fastqc.zip
    screenData=${f%.fastq.gz}_sampled_screen.txt
    echo Found fastq file:$f
    #zcat $f | head -n 400000 > $sampled #make new file with first 100,000 reads
    seqtk sample -s100 $f 100000 > $sampled
    echo "Blasting $sample to look for contaminants"    
    #Lets blast the reads to see what they are...
    head -n 40000 $sampled | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${sampled%stq}
    blastn -db /gpfs/genomes/blast_databases/refseq_rna.00 -query ${sampled%stq} -out ${sampled%_sampled.fastq}_blast.txt -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -max_target_seqs 1 -num_threads 12
    blastn -db /gpfs/genomes/blast_databases/UniVec -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 175000000000 -query ${sampled%stq} -out ${sampled%_sampled.fastq}_vec.txt -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -max_target_seqs 1 -num_threads 12
    awk -F '\t' '{print $8}' ${sampled%_sampled.fastq}_blast.txt | sort | uniq -c | sort -nr | head -n 5 > ${sampled%_sampled.fastq}_hist.txt
    awk -F '\t' '{print $10}' ${sampled%_sampled.fastq}_vec.txt | sort | uniq -c | sort -nr | head -n 5 >> ${sampled%_sampled.fastq}_hist.txt
    echo "Top hits were:"
    cat ${sampled%_sampled.fastq}_hist.txt
    fastqc --noextract $sampled  #run fastqc on first 400,000
    unzip -q -o -d $dir $zip  #unzip fastqc file
    duplicate=$(grep 'Total Deduplicated Percentage' $data) #get Duplicate Percentage line from fastqc_data.txt
    percent=$(echo $duplicate | awk -F'.' '{print $1}' | awk -F' ' '{print $NF}') #get the percentage from that line
    percent=$((100-percent))
    echo "Total Duplicate Percentage for Sample $f [$percent%]" 
    fastq_screen --aligner bowtie --subset 0 --outdir $dir --force --quiet --threads 12 $sampled
    if [ -f "$screenData" ]
    then
          cat $screenData | awk '{print $4}' | while read line #go through the %unmapped column of fastqscreen results
          do
           unmapped=$(grep -Eo '[0-9\.]+')
           for u in $unmapped; do #for each %unmapped, get the %mapped
            mapped="$(echo "100.00-$u" | bc)"
            if [ 1 -eq "$(echo "$mapped >= 5" | bc)" ]; then #if %mapped is greater/equal to 5, then record it
             ref=$(grep -m 1 $u $screenData | awk '{print $1}')
             mapPercents="$mapPercents $ref $mapped%"
            fi
          done
          echo Percent: $percent
          sample=${f%.fastq.gz}
          echo "${sample##*/} : $percent% duplication : $mapPercents" >> $dir/QC/MoreStats.txt #add %duplication and %mapped to the quality report
          done
    else
        echo "$screenData did not exist, so we didn't try to process it to avoid infinite loops and errors" 
    fi
    cat ${sampled%_sampled.fastq}_hist.txt >> $dir/QC/MoreStats.txt
    mv ${dir}/*.html ${dir}/QC
    rm -rf $sampled ${f%.fastq.gz}_sampled_fastqc* ${sampled%stq} ${dir}/*.txt
    echo "Cleaned up temporary files for $f"
  done #for each fastq file
done
fi
