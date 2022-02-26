#!/bin/bash

# automated demultiplexing, QC, backup, cloudbackup, and reporting script.
# Last edited on 4/17/2018

if [[ $# -lt 1 ]]
then
	echo "Usage: startDemultiplexer_v2.sh DirectoryContainingNewlySequencedRunFolders email1 email2 email3 ..."
	echo "Runs continuously so make sure you have screen activated! Errors? Contact Salk IGC"
	here=/gpfs/smb/illumina/runfolders 
	email='nhah@salk.edu mshokhirev@salk.edu ouyang@salk.edu' #if more than one email, do this --> lliu@salk.edu maxshok@gmail.com another@salk.edu another@salk.edu ...

else
	here=$1
	shift 1
	email=$@
fi
outdir=/gpfs/smb/illumina/fastq

echo -e "\nStarting Lucas's IGC Demultiplexer, version 2.3\nScanning $here for new sequencing runs.\nWill demultiplex, run basic QC, and send out email reports to: $email \n "
if ! command -v bcl2fastq &> /dev/null
then
    echo "bcl2fastq could not be found. Please source /gpfs/tools/.bashTools"
    exit
fi
if ! command -v DemultiplexVisualizer &> /dev/null
then
    echo "Demultiplexvisualizer could not be found. Please source /gpfs/tools/.bashTools"
    exit
fi
if ! command -v fastqc &> /dev/null
then
    echo "fastqc could not be found. Please source /gpfs/tools/.bashTools"
    exit
fi
if ! command -v fastqScreen &> /dev/null
then
    echo "fastqScreen could not be found. Please source /gpfs/tools/.bashTools"
    exit
fi
while :
do

####### THIS WHILE LOOP IS ACTUALLY JUST FOR SELECTING THE NEXT AVAILABLE NEW RUNFOLDER! ########
while true; do
	runDirs=$(find $here -maxdepth 2 -name RTAComplete.txt -printf '%h\n') #Find RTAComplete.txt
	#echo Found RTAComplete here: $runDirs
	for dir in $runDirs; do #go through all directories with RTAComplete

	#	echo dir: $dir
		#echo basen: $(basename $dir)
		result=${dir##*/}
	#	echo result: $result
 		#if the directory has a sample sheet AND has not been demultiplexed, then continue with script
		sheetCount=$(find $dir -maxdepth 1 -name "*.csv"| wc -l)
 		if [ $sheetCount -gt 0 ] && [ ! -f $dir/demultiplexed.txt ]; then
  			found=true
  			break
 		#else
  		#	echo RTAComplete/SampleSheet/Unfinished Status not found.
 		fi
	done
	sleep 1
	#echo "Test"
	if [ "$found" = true ]; then #continue with script
 		break
	fi
done

############### AT THIS POINT WE HAVE SELECTED A NEW RUNFOLDER THAT NEEDS DEMULTIPLEXING ###################
echo Unprocessed run folder found: $dir 

startTime=`date +%s`
#Determine whether SE or PE
#sheet=$(find $dir -maxdepth 3 -name *.csv -type f)
sheet=$(find $dir -maxdepth 1 -name "*.csv" -type f)
echo sheet: $sheet
 numN=$(grep -c 'IsIndexedRead="N"' $dir/RunInfo.xml) 
 dateline=$(grep 'Date' $dir/RunInfo.xml)
 date=$(echo $dateline | awk -F"</Date" '{print $1}' | awk -F"Date>" '{print $NF}') #Get date from RunInfo	
 #output = date+pe/se+length
 #output=$date
 output=$(basename $dir | awk -F"_" '{print $1}')
 echo outp: $output
 if [ $numN -ge 2 ]; then #Determine whether SE or PE
  output+="_PE"
 else
  output+="_SE"
 fi

#DETERMINE LENGTH
	if [ $numN -ge 2 ]; then
	 declare -i length=0
	 line=$(grep -m 1 'Read="N"' $dir/RunInfo.xml | tail -n 1 )
	 formatted=$(echo $line | awk -F " />" '{print $1}' | awk -F "<" '{print $NF}')
	 thisLength=$(echo $formatted | awk -F '" IsIndexed' '{print $1}' | awk -F 'Cycles="' '{print $NF}')
	 length=$thisLength
	 output+="${length}"
	else #se
	  line=$(grep -m 1 'NumCycles=' $dir/RunInfo.xml)
	  formatted=$(echo $line | awk -F" />" '{print $1}' | awk -F"<" '{print $NF}')
	  length=$(echo $formatted | awk -F'" IsIndexed' '{print $1}' | awk -F'Cycles="' '{print $NF}')
	  output+="${length}"
	fi
	echo Output folder: $output
	echo Run folder: $dir
	#sheetname=${sheet}
	sheetname=$(basename $sheet)
        echo SampleSheet name: $sheetname
	#cat $sheet
echo "In 3 seconds: Starting bcl2fastq on $dir"
sleep 3
#run bcl2fastq, if no problem then continue with quality check and email, if problem then email error log, add demultiplexed status, and return to beginning
if [ -e ${outdir}/${output} ]
then
	echo "The output directory ${outdir}/${output} already exists! To prevent duplicate reads, we will first delete the entire output folder!"
	for i in $(seq 10 -1 1)
	do
		echo "Deleting ${outdir}/${output} in $i seconds"
		sleep 1
	done 
	rm -r ${outdir}/${output}
	echo "Started demultiplexing $dir and outputing to $outdir/$output using samplesheet $dir/$sheetname after deleting what was already there..."
else
	echo "Started demultiplexing $dir and outputing to $outdir/$output using samplesheet $dir/$sheetname"
fi


read_one=$(grep 'Read Number="1"' $dir/RunInfo.xml | awk -F'"' '{print $4}') #novaseq 10x
echo "Read 1 length:" $read_one

read_two=$(grep 'Read Number="2"' $dir/RunInfo.xml | awk -F'"' '{print $4}') #novaseq 10x
echo "Read 2 length:" $read_two

read_three=$(grep 'Read Number="3"' $dir/RunInfo.xml | awk -F'"' '{print $4}') #novaseq 10x
echo "Read 3 length:" $read_three

echo outdir: $outdir
echo output: $output
echo both: $outdir/$output

if [[ $read_one == "28" && $read_two == "8" && $read_three == "91" ]]; then
        echo "Read lengths indicate this as a 10X run, converting to fastq using 10X options... output is saved here: $dir/log.txt"
	bcl2fastq --use-bases-mask Y28,I8,Y91 --create-fastq-for-index-reads --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8 --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls -r 6 -w 6 -R $dir -o ${outdir}/${output} --interop-dir ${dir}/InterOp --sample-sheet ${dir}/${sheetname} 2> $dir/log.txt
	if [ $? -ne 0 ]
	then	
		echo "ERROR: "
		cat $dir/log.txt
	else
        	echo "Finished demultplexing with 10X options"
	fi

else
        echo "Started demultiplexing using default bcl2fastq parameters... output is saved here: $dir/log.txt"
        bcl2fastq -l DEBUG -r 4 -p 16 -w 4 --barcode-mismatches 0 --no-bgzf-compression --no-lane-splitting -R $dir -o ${outdir}/${output} --sample-sheet ${dir}/${sheetname} 2> $dir/log.txt
	if [ $? -ne 0 ]
	then	
		echo "ERROR: "
		cat $dir/log.txt
	else
        	echo "Finished demultiplexing with default bcl2fastq parameters"
	fi
fi


retval=$?
echo "Demux returned: [$retval]"
if [ $retval -eq 0 ] ; then
 echo Demultiplexing of $dir completed successfully
 touch $dir/demultiplexed.txt #mark run directory as demultiplexed
	#Run graphing programs on DemultiplexingStats.xml
	DemultiplexVisualizer ${outdir}/${output}/Stats/DemultiplexingStats.xml
	echo Graphs created for ${outdir}/${output}
	#Create string that will be used for email attachment at the end
	for i in $(find ${outdir}/${output}/Stats/ -type f -name *.png); do pngs="$pngs -a $i"; done;
	datetime=$(date)
	#get lane.html
        lane=${outdir}/${output}/Reports/html/*/default/Undetermined/unknown/lane.html
	NEWLINE=$'\n'

#Create quality report
 echo More Statistics from FastQC and FastQScreen for "${outdir}/${output}" > ${outdir}/${output}/MoreStats.txt
 echo Blank = Sample did not map at least 5% to any reference genomes >> ${outdir}/${output}/MoreStats.txt
 echo "" >> $outdir/$output/MoreStats.txt

 echo "<HTML><BODY><h1><b>Salk NGS Core</b></h1><nbsp>Sorry this part is not accessible. <nbsp>Please make sure you type in the full link to your sequencing data, or contact the <a href=\"http://www.salk.edu/science/core-facilities/next-generation-sequencing/\">NGS core</a> for help.</BODY></HTML>" > ${outdir}/${output}/index.html

for d in ${outdir}/${output}/*/ ; do #for each project
 base=$(basename $d)

 if [ $base != "Reports" ] && [ $base != "Stats" ] && [ "$base" != *"Undetermined"* ] && [ "$base" != *"QC" ] ; then
   mkdir ${d}QC
   echo Project: $base >> ${outdir}/${output}/MoreStats.txt
   echo Project: $base > ${d}QC/QC.txt
   echo Found project directory: $d

  #static file extensions

  fastqfiles=$(find ${d}*fastq.gz)
  for f in $fastqfiles; do #for each fastq file
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

    fastqc --noextract -o $d $sampled  #run fastqc on first 400,000
    unzip -q $zip -d $d #unzip fastqc file
    duplicate=$(grep 'Total Deduplicated Percentage' $data) #get Duplicate Percentage line from fastqc_data.txt
    percent=$(echo $duplicate | awk -F'.' '{print $1}' | awk -F' ' '{print $NF}') #get the percentage from that line
    percent=$((100-percent))
    echo "Total Duplicate Percentage for Sample $f [$percent%]" 
    fastq_screen --aligner bowtie --subset 0 --outdir $d $sampled
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
	  mv ${f%.fastq.gz}*screen*.html ${f%.fastq.gz}*fastqc.html ${d}QC
	  rm ${f%.fastq.gz}*screen*.txt
	  echo Percent: $percent
	  sample=${f%.fastq.gz}
	  echo "${sample##*/} : $percent% duplication : $mapPercents" >> $outdir/$output/MoreStats.txt #add %duplication and %mapped to the quality report
	  echo "${sample##*/} : $percent% duplication : $mapPercents" >> ${d}QC/QC.txt #add %duplication and %mapped to the quality report
	  done
    else
	echo "$screenData did not exist, so we didn't try to process it to avoid infinite loops and errors" 
    fi
    rm -rf $sampled ${f%.fastq.gz}_sampled_fastqc* ${sampled%stq} 
    echo "Cleaned up temporary files for $f"
  done #for each fastq file
  echo >> $outdir/$output/MoreStats.txt
  echo "Top 5 Organisms and top 5 vectors" > ${d}Top10Organisms.txt
  for i in $(ls ${d}*hist.txt); do printf "${i##*/}\t" >> ${d}Top10Organisms.txt ; done
  printf "\n" >> ${d}Top10Organisms.txt
  paste -d'\t' ${d}*hist.txt  >> ${d}Top10Organisms.txt
  rm ${d}*hist.txt ${d}*_blast.txt ${d}*_vec.txt 
  echo ${d} >> $outdir/$output/MoreStats.txt
  cat ${d}Top10Organisms.txt >> $outdir/$output/MoreStats.txt
  rm ${d}Top10Organisms.txt;
  fi
done
#All of the quality checks for the Undetermined fastq file
undet=$(find ${outdir}/${output}/*R1*fastq.gz)
if [ -f "$undet" ]
then
	name=$(basename $undet)

	#echo "Undetermined: $name" >> $outdir/$output/MoreStats.txt
	 zcat $undet | head -n 400000 > ${undet%.fastq.gz}_sampled.fastq
	 fastqc --noextract -o $outdir/$output ${undet%.fastq.gz}_sampled.fastq
	 unzip -q ${undet%.fastq.gz}_sampled_fastqc.zip -d $outdir/$output

	duplicate=$(grep 'Total Deduplicated Percentage' ${undet%.fastq.gz}_sampled_fastqc/fastqc_data.txt)
	percent=$(echo $duplicate | awk -F'.' '{print $1}' | awk -F' ' '{print $NF}')
	percent=$((100-percent))
	echo "Total Duplicate Percentage for Sample $undet [$percent%]"

	fastq_screen --outdir $outdir/$output --aligner bowtie --subset 0 ${undet%.fastq.gz}_sampled.fastq
	cat ${undet%.fastq.gz}_sampled_screen.txt | awk '{print $4}' | while read line
	do
	 unmapped=$(grep -Eo '[0-9\.]+')
	 for u in $unmapped; do
	 	mapped=$(echo "100.00-$u" | bc)
		comparison=$(echo "$mapped >= 5" | bc)
	 	if [[ 1 -eq "${comparison:-0}" ]]
		 then
	 	    ref=$(grep -m 1 $u ${undet%.fastq.gz}_sampled_screen.txt | awk '{print $1}')
	 	    mapPercents="$mapPercents $ref $mapped%"
	 	    echo "mapPercent: $mapPercents"
	 	fi
	 done
	 echo "$name : $percent% duplication : $mapPercents" >> ${outdir}/${output}/MoreStats.txt
	done
else
	echo "There was no undetermined fastq file found in ${outdir}/${output}"
	echo "There was no undetermined fastq file found in ${outdir}/${output}" >> ${outdir}/{$output}/MoreStats.txt
fi

#checksums
for d in $outdir/$output/*/ ; do
				   base=$(basename $d)
				   if [ $base != "Reports" ] && [ $base != "Stats" ] && [[ "$base" != *"Undetermined"* ]]; then #for each project directory
				   	echo d: $d 
				   	touch ${d}checksum.txt #create checksum file
				   	echo $base fastq checksums >> ${d}checksum.txt #title checksum file
				   	echo  ${NEWLINE} >> ${d}checksum.txt
				   	fastqs=$(find $d*fastq.gz) #get fastq files
				   	for f in $fastqs; do #for each fastqfile
				   		checksum=$(md5sum $f) #checksum the file
				   		echo $checksum >> ${d}checksum.txt #add checksum to checksum file
				   		echo $f complete.
				   	done
				   	echo Checksum for $base done.
					#mkdir ${d}QC
				   fi
done
				 
other=$(find $outdir/$output/*R1*fastq.gz) #do the checksum for Undetermined file
touch $outdir/$output/Undetermined_checksum.txt
checksum=$(md5sum $other)
echo $checksum >> $outdir/$output/Undetermined_checksum.txt
echo $base done.
#mkdir $outdir/$output/QC
mv $outdir/$output/*.txt $outdir/$output/*.html $outdir/$output/QC
rm -r $outdir/$output/*.zip $outdir/$output/*fastqc* 

#echo "Changing permissions on $outdir/$output to 751"
#chmod 751 $outdir/$output

#This is now being handled by a daily Cron Job
#backup to boll 
echo "backing up $outdir/$output to $backup "
rsync -av --progress --ignore-existing $outdir/$output $backup | tail -n 10 > $outdir/$output/backup.log
cat $outdir/$output/backup.log

#automatic cloud backup 
#echo "Starting OracleCloud backup for $output"
#tar czf ${outdir}/${output}.tar.gz $outdir/$output 
#echo "Uploading $outdir/${output}.tar.gz to igc-archive"
#java -jar /dev/gpfs/tools/OracleCloud/ftmcli.jar upload -a -I a429964 -auth-url https\://us2.storage.oraclecloud.com -service Storage -user igc@salk.edu igc-archive $outdir/${output}.tar.gz | tail -n 5 > $outdir/${output}/cloudbackup.log
#backupToOracleCloud $outdir | tail -n 10 > $outdir/$output/cloudbackup.log
#cloudgood=$(grep successful $outdir/${output}/cloudbackup.log | wc -l )
#echo "Removing local copy of $outdir/${output}.tar.gz"
#rm $outdir/${output}.tar.gz

#Email all of it to specified emails
tail -n 1000 $dir/log.txt > $dir/taillog.txt
endTime=`date +%s`
runtime=$(( (endTime-startTime)/60 ))
unset links
for i in $outdir/$output/*_*/ 
do
	links="${NEWLINE}igc2.salk.edu/illumina/runs/$output/$(basename $i) ${links}"
done
#if [ $cloudgood -gt 0 ]
#then
	mail -a $dir/taillog.txt -a $outdir/$output/MoreStats.txt -a $lane $pngs -s "$output Complete" $email <<< "Completed on ${datetime}. Took $runtime mins.  ${NEWLINE}Statistics from FastQC, FastQScreen, and blast attached in MoreStats.txt ${NEWLINE}Flowcell/Lane Summaries and Top Unknown Barcodes attached in lane.html. ${NEWLINE} ${NEWLINE} Project links : $links ${NEWLINE} ${NEWLINE}"
#else 
#	mail -a $dir/taillog.txt -a $outdir/$output/MoreStats.txt -a $lane $pngs $outdir/${output}/cloudbackup.log -s "$output Complete" $email <<< "Completed on ${datetime}. Took $runtime mins.  ${NEWLINE}Statistics from FastQC, FastQScreen, and blast attached in MoreStats.txt ${NEWLINE}Flowcell/Lane Summaries and Top Unknown Barcodes attached in lane.html. ${NEWLINE} ${NEWLINE} Project links : $links ${NEWLINE} ${NEWLINE}"
#fi

echo Emails sent.
touch $dir/demultiplexed.txt #mark demultiplexed status
unset pngs
unset found
else #if demux failed, then email error log to specified emails
 tail -n 100 $dir/log.txt > $dir/taillog.txt
 cat $dir/log.txt
 mail -a $dir/taillog.txt -s "$output Demux Failed" $email <<< "Demultiplexing for $output failed." #email error log
 touch $dir/demultiplexed.txt #mark demultiplexed status
 echo Demultiplexing error encountered, restarting script.
 unset found
fi

done
