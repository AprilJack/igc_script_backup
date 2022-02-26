/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * (M)RNA (A)SSEMBLY for (P)ROTEOGENOMIC(S), MAPS, is a general transcript java assembly tool which takes aligned reads and outputs a diverse library of transcripts. 
 * Unlike extant transcript assembly tools which provide parsimonious assembly, MAPS tries to maximize the diversity of the transcript library generated in
 * order to facilitate downstream searching for rare or un-annotated transcripts. While MAPS assembles transcripts, it outputs peptide-level assemblies as well. 
 * MAPS outputs a GTF file containing all of the 3-frame translated ORFs within a given size range (For peptide searching). It also outputs basic abundance and depth stats 
 * for each ORF/transcript. In addition, if -consensus is specified, MAPS will build a "consensus" sequence from the reads, enabling handling of sample-specific point mutations.
 * MAPS handles stranded and unstranded datasets (guesses strand based on longest ORF) and works with paired or unpaired libraries. MAPS allows the diversity of the libraries to be 
 * adjusted by varying the stringency parameter. Finally, MAPS can be run in parallel to speed up assembly, or in low memory mode (-p 1) if RAM is limited to < 16 GB.
 * 
 * MAPS was designed for detecting short ORFs in RNA-Seq data at the Salk Institute.
 *  
 * @author Maxim Nikolaievich Shokhirev (C) 2015-2016
 */
public class MAPS implements Comparator<int[]>{

	double version = 1.6;
	HashMap<String,LinkedList<int[]>> reads = new HashMap<String,LinkedList<int[]>>();
	HashMap<String,LinkedList<Transcript>> transcripts = new HashMap<String,LinkedList<Transcript>>();
	HashMap<String,LinkedList<Exon>> exons = new HashMap<String,LinkedList<Exon>>();  						//forward strand exons
	HashMap<String,HashMap<Integer,Integer>> junctions = new HashMap<String,HashMap<Integer,Integer>>(); 	//exon junctions (stranded)
	
	HashMap<String,String> nucToPep = buildNucToPep(); 
	
	int totalFragments = 0;	
	double avgFragLength=0;
	int totalSNPs = 0;
	int ext = 100;   							//extend all reads by this number of bps in each direction
	HashMap<String,Integer> totalORFs = new HashMap<String,Integer>();
	HashMap<String,HashMap<String,int[]>> readMap = new HashMap<String,HashMap<String,int[]>>();   //need to save read name for paired!
	double mult = 1; 							//multiplier used to figure out coverage thresholds during filtering
	HashMap<String,byte[]> chrs = null;
	HashMap<String,Integer> chrLengths = new HashMap<String,Integer>();
	boolean verbose = false;  					//set if you want to see a bunch of debug output
	boolean filter = false;						//set if you want to filter exons after agglomeration by minimum length and read count
	int cores = Runtime.getRuntime().availableProcessors()-1;
	String prefix = "";							//output with this prefix
	int fmax = 1000;  							// max length of a transcript to report as "short" (saved in a separate transcript file)
	int fmin = 8;  								// max length of a transcript to report as "short" (saved in a separate transcript file)
	long avgLength = 0;
	ArgParser ap = null;
	int libType = 0;  //fr, rf, ff, rr
	boolean unstranded = false;
	boolean paired = false;
	double stringency = 1;
	double overlapMax = 0.0;
	boolean buildConsensus = false;
	boolean keepTempFiles = false;
	boolean extendTranscripts = true;
	int minMapQ = 0;
	
	HashMap<String,LinkedList<SNP>> snps = new HashMap<String,LinkedList<SNP>>();
	
	Pattern p = Pattern.compile("([0-9]+[MIDNSHP=X])");
	
	public MAPS(String[] args) {
		ap = new ArgParser("MAPS v"+version+" uses aligned reads from unsorted SAMs to assemble a diverse transcriptome, and outputs all 3frame ORFs and stats."
				+ "\nHas options for controlling library diversity, using the read sequence to account for mutations, and basic coverage/FPKM/TPM calculation."
				+ "\nDeveloped by Max Shokhirev at the Salk Institute IGC Core (C) 2015.\n"
				+ "\nUsage: MAPS alignment.sam [alignment2.sam] ... [alignmentN.sam] [opts]\n"
				+ "\n\nAll exon-exon junctions, exons, SNPs, args, 3frame translated fasta, nucleotide fasta, and annotations are printed to prefix.junc/exon/SNPs/args/pep/nuc/gtf files, respectively."
				+ "\nThe sequence file used during mapping is required for ORF-level analyses, extension of transcripts up/down stream, consensus sequence construction, and low-memory mode.");
		ap.registerArg("seq",  null, "A fasta file containing the genomic sequences of the chromosomes that the reads were mapped to. Keep chr names consistent.");
		Date date = new Date();
        SimpleDateFormat format =
                new SimpleDateFormat("dd-MMM-yyyy;HH-mm");
		ap.registerArg("o", format.format(date),"Output prefix used for outputing transcript files");
		ap.registerArg("p", ""+(Runtime.getRuntime().availableProcessors()/2),"The number of cores to use for read agglomeration. Run with -p 0 to turn on low-memory mode (high disk usage and slower).");
		ap.registerArg("library","fr", "Uses library design to deconvolute strand identity of reads (un,f,r,fr,rf,ff,rr)");		
		ap.registerArg("stringency","0.5","Value used to control the diversity of the assembled transcriptome (decrease to increase transcriptome diversity/size). Using a very low value can result in combinatorial explosion! Use 0 to cycle through 0.1 to 1.0 in steps of 0.9");
		ap.registerArg("consensus","false","Uses a read voting scheme that adjusts for depth to account for mutations. Generates SNP file in pgSNP format. Adds notation to ORFs.");
		ap.registerArg("v", "false", "Verbose. Set if you would like to see additional messages/stats during the run.");
		ap.registerArg("fmax","1000000","ORFs will be generated with lengths <= fmax peptides");
		ap.registerArg("fmin","10","ORFs will be generated with lengths => fmin peptides");
		ap.registerArg("minMAPQ","0","Reads with MAPQ less than this will not be used!");
		ap.registerArg("keepTempFiles","false","If running in low-memory mode (by setting -p 0), set to true to keep the temporary chr-separated fasta and sam files.");
		ap.registerArg("printExons","true","Outputs assembled exons as a gff file.");
		ap.registerArg("printJunctions","true","Outputs junctions as a gff file.");
		ap.registerArg("mult", "0", "Coverage depth multiplier. Leave at 0 to automatically adjust by totalReads/10^6. Increase if skipping chromosomes.");
		ap.registerArg("ext", null, "The amount each read will be extended. Set to override the default of 300(1.0-stringency^2)");
		ap.registerArg("noTranscriptExtension","false","MAPS will automatically extend transcripts until a stop codon is found in all frames and both directions.");
	
		ap.parseArgs(args);
		if(ap.getList().size() == 0)
			ap.printUsage();
		long tic = System.currentTimeMillis();
		fmax= ap.getAsInt("fmax");
		fmin =ap.getAsInt("fmin");
		
		if(ap.get("library").compareTo("fr")==0)
		{
			libType = 0;
			paired = true;
		}
		else if(ap.get("library").compareTo("rf")==0)
		{
			libType = 1; paired = true;
		}
		else if(ap.get("library").compareTo("ff")==0)
		{
			libType = 2; paired =true;
		}
		else if(ap.get("library").compareTo("rr")==0)
		{
			libType = 3; paired =true;
		}
		else if(ap.get("library").compareTo("un")==0)
			unstranded=true;
		else if(ap.get("library").compareTo("r")==0)
			libType = 1;   
		verbose = ap.getAsBoolean("v");
		
		cores = Math.min(cores,ap.getAsInt("p"));
		prefix=ap.get("o");
		buildConsensus = ap.getAsBoolean("consensus");
		filter = true; //ap.getAsBoolean("filter");
		stringency = Math.max(0,Math.min(1,ap.getAsDouble("stringency")));
		ext = (int) (300*(1.0-Math.pow(stringency,2)));
		extendTranscripts = !ap.getAsBoolean("noTranscriptExtension");
		
		overlapMax = 0.00; //Math.max(0,Math.min(1,ap.getAsDouble("overlap")));
		keepTempFiles=ap.getAsBoolean("keepTempFiles");
		PrintWriter argWriter = null;
		try{
			argWriter = new PrintWriter(prefix+".args.txt");
			{
				String line = "";
				for(String arg: args)
					line+=" "+arg;
				argWriter.println("MAPS"+line);
				for(String arg:ap.getAllArgs())
				{
					argWriter.println(arg);
				}
			}
		}catch(Exception e)
		{e.printStackTrace();}
		finally{
			argWriter.close();
		}
		PrintWriter gtf = null;
		PrintWriter full_gtf = null;
		try {
			 if(ap.get("seq")!=null)
				 gtf= new PrintWriter(new GZIPOutputStream(new FileOutputStream(prefix+"_ORFs.gtf.gz")));
			 full_gtf= new PrintWriter(new GZIPOutputStream(new FileOutputStream(prefix+".gtf.gz")));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		if(ap.get("seq") == null || !(new File(ap.get("seq")).exists()))
		{
			System.err.println("Did not detect proper fasta sequence file. Please specify with -seq if you want ORF-level analysis, consensus SNP calling, and/or low-memory mode operation.");
			cores = Math.max(cores, 1);
			buildConsensus = false;
		}
		if(ap.get("ext") != null)
		{
			System.err.println(String.format("Extending each read by %d (default for stringency = %3.2f is %d)",ap.getAsInt("ext"),stringency,ext));
			ext = ap.getAsInt("ext");
		}
		
		if(ap.get("minMAPQ")!= null)
		{
			minMapQ = ap.getAsInt("minMAPQ");
		}
		System.err.println((new Date()).toString()+":Starting MAPS. Stranded="+!unstranded+" Paired="+paired+" LibType="+ap.get("library")+" Consensus="+buildConsensus+" Extending reads: "+ext+" extendTranscripts: "+extendTranscripts+". Stringency = "+stringency+" on "+cores+" threads. Min MAPQ: "+minMapQ+". Output will go to "+prefix);
		long tic2 = System.currentTimeMillis();
		if(cores == 0)
		{
			cores =1;
			//lets read in the sam file and save it as a bunch of temp chr sam files
			//also lets read in the seq file and break it up by the chrs
			LinkedList<String> chrNames = splitSeq(ap.get("seq"));
			splitSams(ap.getList(),chrNames);
			for(String chrName: chrNames)
			{
				if(verbose)System.err.println("Working with "+chrName);
				tic2 = System.currentTimeMillis();
				readInChr(chrName,prefix+"_"+chrName+".fa");
				if(verbose) System.err.println("Loading genome sequence took "+String.format("%3.3f",((System.currentTimeMillis()-tic2)/60000.0))+" minutes.");
				LinkedList<String> theSamList= new LinkedList<String>();
				theSamList.add(prefix+"_"+chrName+".sam");
				try {
					loadReads(theSamList);
				} catch (NumberFormatException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
				if(verbose) System.err.println("Loading reads took "+String.format("%3.3f",((System.currentTimeMillis()-tic2)/60000.0))+" minutes.");
				summarizeReads();
				sortReads();
				
				//if(verbose) System.err.println("Building exons/junctions and filtering took "+String.format("%3.3f",(System.currentTimeMillis()-tic2)/60000.0)+" minutes.");
				if(verbose && !keepTempFiles)
				{
					File fa = new File(prefix+"_"+chrName+".fa");
					fa.delete();
					File samFile = new File(prefix+"_"+chrName+".sam");
					samFile.delete();
					System.err.println("Deleted temporary sam and sequence file for "+chrName);
				}
			}  
			
		}
		else
		{
			if(ap.get("seq") != null && (new File(ap.get("seq")).exists()))
				readInAllChrs(ap.get("seq"));
			//if(verbose) System.err.println("Loading genome sequence took "+String.format("%3.3f",((System.currentTimeMillis()-tic2)/60000.0))+" minutes.");
			tic2 = System.currentTimeMillis();
			try {
				loadReads(ap.getList());
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (IOException e) {

				e.printStackTrace();
			}
			if(verbose) System.err.println("Loading reads took "+String.format("%3.3f",((System.currentTimeMillis()-tic2)/60000.0))+" minutes.");	
			System.err.println();
			System.err.println("MAPS starting read agglomeration into transcripts using "+cores+" threads...");
			summarizeReads();  	//output stats on the reads
			tic2 = System.currentTimeMillis();
			sortReads();		//sort the reads by their genomic coords
			//if(verbose) System.err.println("Sorting the reads took "+String.format("%3.3f",((System.currentTimeMillis()-tic2)/60000.0))+" minutes.");
		}		//reads->exons+junctions
		
		processReads();
		if(verbose) System.err.println("Building exons/junctions and filtering took "+String.format("%3.3f",(System.currentTimeMillis()-tic2)/60000.0)+" minutes.");

		reads = null;
		try {
			if(ap.getAsBoolean("printJunctions"))
				outputJunctions();	//output all junctions to a gtf file for debuging purposes (can be loaded into genome browser)
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			if(ap.getAsBoolean("printExons"))
				outputExons();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		PrintWriter pw = null;
		PrintWriter nucFasta = null;
		if(ap.get("seq")!= null)
		{
			try {
				pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(prefix+".pep.gz")));
				nucFasta = new PrintWriter(new GZIPOutputStream(new FileOutputStream(prefix+".nuc.gz")));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		int totalTranscripts = 0;
		HashSet<String> chrSet = new HashSet<String>();
		for(String chrStr: exons.keySet())
		{
			if(exons.get(chrStr).size() > 10)
			{
				chrSet.add(chrStr.substring(0,chrStr.length()-1));
				totalTranscripts+=transcripts.get(chrStr).size();
			}
		}
		//if(verbose)System.err.println("Outputing transcripts...");
		tic2 = System.currentTimeMillis();
		if(!verbose) System.err.println("Outputing transcripts:");
		if(!verbose) System.err.println("|0%                     25%                      50%                      75%                  100%|");
		int transcriptsSoFar = 0;
		for(String chr: chrSet)
		{ 
			if(transcripts.get(chr+"+")!= null)
			{
				outputTranscripts(transcripts.get(chr+"+"),chr+"+",true,pw,gtf,full_gtf,nucFasta);
				transcriptsSoFar+=transcripts.get(chr+"+").size();
			}
			if(!verbose)
			{
				String bar = "";
				int percentDone = (int)((100.0*transcriptsSoFar)/totalTranscripts);
				for(int i = 0; i < percentDone; i++)
					bar+="=";
				System.err.print(bar+"                                               \r");
			}
			if(!unstranded) //stranded library
			{
				if(transcripts.get(chr+"-")!= null)
				{
					outputTranscripts(transcripts.get(chr+"-"),chr+"-",false,pw,gtf,full_gtf,nucFasta);
					transcriptsSoFar+=transcripts.get(chr+"-").size();
				}
				if(!verbose)
				{
					String bar = "";
					int percentDone = (int)((100.0*transcriptsSoFar)/totalTranscripts);
					for(int i = 0; i < percentDone; i++)
						bar+="=";
					System.err.print(bar+"                                               \r");
				}
			}

			
		}
		full_gtf.close();
		if(ap.get("seq")!= null)
		{
			pw.close();
			gtf.close();
			nucFasta.close();
		}
		int totalOverall = 0;
		if(ap.get("seq")!= null)
		{
			for(String chr: totalORFs.keySet())
			{
				if(verbose)System.err.println(chr+" had "+totalORFs.get(chr)+" ORFs");
				totalOverall+=totalORFs.get(chr);
			}
		}
		double mins = ((int)((System.currentTimeMillis()-tic)/6000.0))/10.0;
		if(!verbose)System.err.println();
		if(ap.get("seq")!= null)
			System.err.println("Finished! Took: "+mins+" mins to process: "+(totalFragments/1000000)+" million fragments resulting in "+totalTranscripts+" transcripts with "+totalOverall+" ORFs ("+fmin+"-"+fmax+" aa). Output prefix: "+prefix);
		else
			System.err.println("Finished! Took: "+mins+" mins to process: "+(totalFragments/1000000)+" million Fragments resulting in "+totalTranscripts+" transcripts. Output prefix: "+prefix);
	}



	private void splitSams(LinkedList<String> sams, LinkedList<String> chrNames)
	{
		long tic = System.currentTimeMillis();
		HashMap<String,StringBuilder> builders = new HashMap<String,StringBuilder>();
		HashMap<String,PrintWriter> printers = new HashMap<String,PrintWriter>();
		try{
			for(String chr: chrNames)
			{
				builders.put(chr,new StringBuilder(1000000));
				printers.put(chr, new PrintWriter(prefix+"_"+chr+".sam"));
			}
			int thusFar = 0;
			for(String sam: sams)
			{
				SuperScanner fs = new SuperScanner(sam);
				while(fs.hasMore())
				{
					//lets scan through this sam file and throw the reads around
					String line = fs.getLine();
					if(line.startsWith("@"))
					{
						//header we can ignore...
					}
					else
					{
						//lets save the read for later processing
						String[] split = line.split("\t");
						if(split.length > 9)
						{
							String chr = split[2];
							if(builders.get(chr) != null)
							{
								if(builders.get(chr).length() > builders.get(chr).capacity()*0.9)
								{
									printers.get(chr).println(builders.get(chr).toString());
									builders.put(chr, new StringBuilder(1000000));
								}
								builders.get(chr).append(line);
								builders.get(chr).append(System.lineSeparator());
							}
							thusFar++;
							if(verbose && thusFar%1000000 == 0) System.err.println("Split up "+thusFar+" reads by chr");
						}
					}
				}
				
			}
			for(String chr: chrNames)
			{
				printers.get(chr).println(builders.get(chr).toString());
				printers.get(chr).close();
			}
			if(verbose)System.err.println("Splitting mapped reads took "+((System.currentTimeMillis()-tic)/60000)+" minutes.");
		}catch(Exception e){e.printStackTrace();}
		
	}


	private LinkedList<String> splitSeq(String seqFile)
	{
		LinkedList<String> result = new LinkedList<String>();
		SuperScanner fs = new SuperScanner(seqFile);
		StringBuilder sb = new StringBuilder(30000000);
		String chr = null;
		while(fs.hasMore())
		{
			String line = fs.getLine();
			if(line.startsWith(">"))
			{
				if(chr != null)
				{
					//lets save the chr to a separate file
					try{
						PrintWriter pw = new PrintWriter(prefix+"_"+chr+".fa");
						pw.println(">"+chr);
						pw.println(sb.toString());
						pw.close();
						result.add(chr);
						if(verbose) System.err.println("Saved temporary fastq file: "+prefix+"_"+chr+".fa");
					}catch(Exception e){} 
				}
				chr = line.substring(1);
				sb = new StringBuilder(30000000);
			}
			else
			{
				sb.append(line); sb.append(System.lineSeparator());
			}
		}
		if(sb != null)
		{
			try{
				PrintWriter pw = new PrintWriter(prefix+"_"+chr+".fa");
				pw.println(">"+chr);
				pw.println(sb.toString());
				pw.close();
				result.add(chr);
				if(verbose) System.err.println("Saved temporary fastq file: "+prefix+"_"+chr+".fa");
			}catch(Exception e){}
		}
		return result;
	}
	
	HashMap<String,short[]> readSeqs = new HashMap<String,short[]>();
	PrintWriter pw = null;
	PrintWriter pw2 = null;
	private void loadReads(LinkedList<String> sams) throws NumberFormatException, IOException {
		if(readMap == null)
		{
			readMap = new HashMap<String,HashMap<String,int[]>>();   //need to save read name for paired!
		}
		if(buildConsensus && (pw == null || pw2 == null))
		{

			try {
				pw = new PrintWriter(prefix+".SNPs.txt");
				pw2 = new PrintWriter(prefix+".allSNPs.txt");
				pw.println("track type=pgSnp visibility=1 name=\""+prefix+" SNPs\" description=\""+prefix+" SNPs\"");
				pw2.println("PeakID\tchr\tstart\tstop\tstrand\tMutation\ttotalCount\tmutantCount\tmutantFreq");
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		for(String sam: sams)
		{
			if(verbose) System.err.println("Loading sam file: "+sam);
			BufferedReader s = null;
			
			InputStream fileStream = new FileInputStream(sam);
			if(sam.endsWith(".gz"))
				s = new BufferedReader(new InputStreamReader(new GZIPInputStream(fileStream)));
			else
				s = new BufferedReader(new InputStreamReader(fileStream));
			while(s.ready())
			{
				String line = s.readLine();
				if(line.startsWith("@"))
				{
					//header we can ignore...
				}
				else
				{
					//lets save the read for later processing
					String[] split = line.split("\t");
					//System.err.println("L:"+line);
					if(split.length > 9)
					{
						int flag = Integer.parseInt(split[1]);
						String chr = split[2];
						int start = Integer.parseInt(split[3]);
						int mapq = Integer.parseInt(split[4]);
						if(mapq>= minMapQ)
						{
							String cigar = split[5];
							//now lets tally the bases in the seq
							boolean mateUnmapped = false;
							boolean plus = true;
							boolean mateplus = true;
							boolean first = false;
							if(((flag & (1 << 3)) != 0))
								mateUnmapped = true;
							if(((flag & (1 << 4)) != 0))
								plus = false;
							if(((flag & (1 << 5)) != 0))
								mateplus = false;						
							if(((flag & (1 << 6)) != 0))
								first = true;						
							if(paired && (libType==0 || libType == 2))  //fr, ff
							{
								if(!first)
									plus = mateplus;   //assume the first segment is accurate
							}
							else if(paired && (libType==1|| libType==3))  //rf, rr
							{
								if(first)
									plus = !plus;
								else
									plus = !mateplus;		
							}
							else if(libType == 1) //r
							{
								//just need to flip the sign since we are sequencing the complement
								plus = !plus;
							}
							//all other combinations don't care: un,f
							String plusString = "+";
							if(!plus&&!unstranded) plusString = "-";						
							//calculate the start/stop pairs for the read based on cigar
							int[] readRegions = getReadRegions(cigar,start);
							//A=65 C=67 G=71 T=84 (other values are ignored)
							if(buildConsensus && readSeqs.get(chr) != null && split[9].length() > 1)  //building consensus and there is a sequence associated with the read
							{
								int counter = 0;
								String editedSeq = getEditedReadSequence(cigar, split[9]);
								byte[] bytes = editedSeq.getBytes();
								for(int r = 0; r < readRegions.length; r+=2)
								{
									
									for(int i = readRegions[r]; i < readRegions[r+1]; i++ )
									{
										try{
											if(readSeqs.get(chr).length > i*4)
											{
												if(bytes[counter]==65||bytes[counter]==97  && (readSeqs.get(chr)[i*4] & 0xFFFF) < 65535)
													readSeqs.get(chr)[i*4]++;
												else if(bytes[counter]==67||bytes[counter]==99  && (readSeqs.get(chr)[i*4+1] & 0xFFFF) < 65535)
													readSeqs.get(chr)[i*4+1]++;
												else if(bytes[counter]==71||bytes[counter]==103  && (readSeqs.get(chr)[i*4+2] & 0xFFFF) < 65535)
													readSeqs.get(chr)[i*4+2]++;
												else if(bytes[counter]==84||bytes[counter]==116  && (readSeqs.get(chr)[i*4+3] & 0xFFFF) < 65535)
													readSeqs.get(chr)[i*4+3]++;
											}
										}
										catch(Exception e){System.err.println(String.format("ERROR Consensus bytesl=%d counter=%d readSeqsl=%d ix4=%d cigar=%s",bytes.length,counter,readSeqs.get(chr).length,i*4,cigar));}
										counter++;
									}
								}
							}
							if(readMap.get(chr)==null)
							{
								readMap.put(chr,new HashMap<String,int[]>());
							}
							if(readMap.get(chr).get(split[0])== null)  //first time seeing a read with this name
							{
								if(paired && !mateUnmapped)
									readMap.get(chr).put(split[0], readRegions);  //just the first part
								else
								{
	
									//this is it (just one read)
									if(reads.get(chr+plusString) == null) reads.put(chr+plusString, new LinkedList<int[]>());								
									
									avgFragLength+=getFragLength(readRegions);
									reads.get(chr+plusString).add(readRegions);
									readMap.get(chr).remove(split[0]);
									totalFragments++;
									if(totalFragments%100000 == 0)
									{
										System.err.print(String.format("Read in %3.1f million aligned reads. Total memory used: %d MB                 \r",(totalFragments/1000000.0),getGigsUsed()));
									}
								}
							}
							else
							{
								//lets combine the regions and resave them
								int[] previousRegions = readMap.get(chr).get(split[0]);
								LinkedList<int[]> combined = combineReadRegions(previousRegions, readRegions);
								if(reads.get(chr+plusString) == null) reads.put(chr+plusString, new LinkedList<int[]>());
								
								avgFragLength+=getFragLength(combined);
								reads.get(chr+plusString).addAll(combined);
								readMap.get(chr).remove(split[0]);// no longer need the old one so lets recycle it!
								totalFragments++;
								if(totalFragments%100000 == 0)
								{
									System.err.print(String.format("Read in %3.1f million aligned fragments. Total memory used: %d MB                 \r",(totalFragments/1000000.0),getGigsUsed()));
								}
							}
						}//didn't pass mapq check
					}//bad sam file line
				}//non-header
			}//while Scanner ready
			s.close();
		}
		readMap = null;
		if(verbose) System.err.println();	
		else System.err.println("Building consensus sequence");
		//lets create the consensus
		if(buildConsensus)
		{
			byte[] letters = {'A','C','G','T'};
			int peakID = 1;
			for(String chr: readSeqs.keySet())
			{
				if(verbose) System.err.println("Building consensus sequence from reads on "+chr);
				for(int pos = 1; pos < readSeqs.get(chr).length/4; pos++)
				{
					int max = 0;
					int index = 0;
					int index_best_mutant = 0;
					int best_mutant_count = 0;
					int total = 0;

					byte old = chrs.get(chr)[pos-1];
					for(int i = pos*4; i <(pos+1)*4 ; i++) {

						if( (readSeqs.get(chr)[i] & 0xFFFF) > max)
						{
							max = (readSeqs.get(chr)[i] & 0xFFFF);
							index = i;
						}
						if( (readSeqs.get(chr)[i] & 0xFFFF) > best_mutant_count && letters[i%4]!=old)
						{
							index_best_mutant = i;
							best_mutant_count = (readSeqs.get(chr)[i] & 0xFFFF) ;
						}
						//if(letters[i%4] !=old)
						total+=(readSeqs.get(chr)[i] & 0xFFFF);
					}
					if(best_mutant_count > 5 && (100.0*best_mutant_count)/total >= 25)
					{
						//output
						pw2.println(prefix+peakID+"\t"+chr+"\t"+(pos-1)+"\t"+(pos)+"\t+\t"+((char)old)+"/"+((char)letters[index_best_mutant%4])+"\t"+total+"\t"+best_mutant_count+"\t"+((int)((100.0*best_mutant_count)/total)));
						peakID++;
					}
					byte newChar = 65;  //A
					if(index%4 == 1) newChar = 67; //C
					else if(index%4 == 2) newChar = 71;  //G
					else if(index%4 == 3) newChar = 84;  //T
					
					if(max > 4 && newChar !=old)
					{
						double score = 0.5+0.5*Math.exp(-max*0.01);
						if(max/(total*1.0) >= score) 
						{						
							chrs.get(chr)[pos-1]=newChar;
							totalSNPs++;
							pw.println(chr+"\t"+(pos-1)+"\t"+(pos)+"\t"+((char)newChar)+"\t1\t"+total+"\t"+max);
							if(snps.get(chr) == null) snps.put(chr,new LinkedList<SNP>());
							snps.get(chr).add(new SNP(pos-1,(char)old,(char)newChar));
						}
					}
				}
			}
		}
		System.err.println();
		if(pw != null)
			pw.close();
		if(pw2 != null)
			pw2.close();
		if(verbose && !buildConsensus)System.err.println("Done loading "+totalFragments);
		else if(verbose && buildConsensus) System.err.println("Done loading  "+totalFragments+". Detected "+totalSNPs+" SNPs");
		mult = totalFragments/1000000.0;
		double tempMult = ap.getAsDouble("mult");
		if(tempMult != 0)
			mult=tempMult;
		avgFragLength/=totalFragments;
		if(verbose) System.err.println(String.format("Coverage multiplier was set at %3.1f average fragment length was %3.1f",mult,avgFragLength));
	}
	
	private int getFragLength(int[] regions)
	{
		int sum = 0;
		for(int i = 0; i < regions.length; i+=2)
		{
			sum+=(regions[i+1]-regions[i]);
		}
		return sum;
	}
	private int getFragLength(LinkedList<int[]> regions)
	{
		int sum = 0;
		for(int[] region: regions)
		{
			for(int i = 0; i < region.length; i+=2)
			{
				sum+=(region[i+1]-region[i]);
			}
		}
		return sum;
	}

/*
	private String splitSeqByCigar(String string, String cigar) {
		Matcher m = p.matcher(cigar);
		//int puffs = 0;
		String output = "";
		int e = 0;
		while(m.find())
		{
			String g = m.group(1);
			if(g.endsWith("D"))
			{
				int skip =Integer.parseInt(g.substring(0,g.length()-1));
				for(int i = 0; i < skip; i++)
					output+="-";   //output gets some dashes to indicate that we are skipping over some bases in the ref				
				//e+=skip;  //do not record these since they don't have an assembly counterpart
			}
			else if(g.endsWith("I"))
			{
				int skip =Integer.parseInt(g.substring(0,g.length()-1));
				e+=skip;
			}
			else if(g.endsWith("N"))
			{
				int skip =Integer.parseInt(g.substring(0,g.length()-1));
				for(int i = 0; i < skip; i++)
					output+="-";   //output gets some dashes to indicate that we are skipping over some bases in the ref
			}
			else if(g.endsWith("S")||g.endsWith("H"))  //assume that these only occur at the beginning or end of a read...
			{

				//do nothing (according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/figure/F1/, SAM positions ignore these)
				
			}
			else
			{
				int record= Integer.parseInt(g.substring(0,g.length()-1));
				output+=string.substring(e,e+record);
				e+=record;
			}
	//		puffs++;
		}
		return output;
	}
*/

	private void summarizeReads()
	{
		HashMap<String,Integer> readCounts = new HashMap<String,Integer>();
		for(String chr: reads.keySet())
		{
			if(readCounts.get(chr.substring(0,chr.length()-1)) == null) readCounts.put(chr.substring(0,chr.length()-1),reads.get(chr).size());
			else readCounts.put(chr.substring(0,chr.length()-1), readCounts.get(chr.substring(0,chr.length()-1))+reads.get(chr).size());
		}
		for(String chr: readCounts.keySet())
		{
			if(readCounts.get(chr) < 1000)
			{
				reads.remove(chr+"+");
				reads.remove(chr+"-");
				if(readCounts.get(chr)>0 )
					System.err.println(chr+" had only "+readCounts.get(chr)+" read fragments and is ommited from further analysis...");
			}
			else
			{
				if(verbose)System.err.println(chr+" had "+readCounts.get(chr)+" read regions.");
			}

		}
		
	}
	private void readInAllChrs(String filename)
	{
		chrs = new HashMap<String,byte[]>();
		SuperScanner fs = new SuperScanner(filename); 
		String chr = null;
		StringBuilder sb = new StringBuilder(1000000);
		while(fs.hasMore())
		{
			String line = fs.getLine();
			if(line.startsWith(">"))
			{
				if(chr != null)
				{
			
					chrs.put(chr,sb.toString().toUpperCase().getBytes());
					chrLengths.put(chr, chrs.get(chr).length);
					if(buildConsensus)
					{
						readSeqs.put(chr,new short[chrs.get(chr).length*4]);
					}
					System.err.print("Loaded "+chr+" consisting of "+chrs.get(chr).length+" nts. Mem total in MB: "+getGigsUsed()+"             \r");
						chr = line.substring(1);
						sb = new StringBuilder(10000000);
				 }
				 else
				 {
					chr = line.substring(1);
				 }
			 }
			 else
			 {
				 sb.append(line);
			 }
		 }
		 if(chr != null)
		 {
			chrs.put(chr,sb.toString().toUpperCase().getBytes());
			chrLengths.put(chr,chrs.get(chr).length);
			if(buildConsensus)
			{
				readSeqs.put(chr,new short[chrs.get(chr).length*4]);
			}
			System.err.print("Loaded "+chr+" consisting of "+chrs.get(chr).length+" nts. Mem total in MB: "+getGigsUsed()+"             \r");
		 }
		 if(verbose)System.err.println();
		 fs.close();		
	}
	
	private void readInChr(String chrToLoad, String filename)
	{
		 chrs = new HashMap<String,byte[]>();
		 SuperScanner fs = new SuperScanner(filename); 
		 System.err.println("Reading in "+chrToLoad);
		 StringBuilder sb = new StringBuilder(10000000);
		 String chr = null;
		 boolean reading = false;
		 while(fs.hasMore())
		 {
			 String line = fs.getLine();
			 if(line.startsWith(">"))
			 {
				 if(chr != null && chrToLoad.compareTo(chr) ==0)
				 {
					chrs.put(chr,sb.toString().toUpperCase().getBytes());
					chrLengths.put(chr,chrs.get(chr).length);
					if(buildConsensus)
					{
						readSeqs.put(chr,new short[chrs.get(chr).length*4]);
					}
					if(verbose)System.err.print("Loaded "+chr+" consisting of "+chrs.get(chr).length+" nts. Mem total in MB: "+getGigsUsed()+"             \r");
					chr = line.substring(1);
					break;
				 }
				 else
				 {
					 chr = line.substring(1);
					 if(chr.compareTo(chrToLoad) == 0)
						 reading = true;
				 }
			 }
			 else if(reading)
				 sb.append(line);
		 }
		 if(chr != null && chr.compareTo(chrToLoad) ==0)
		 {
			chrs.put(chr,sb.toString().toUpperCase().getBytes());
			chrLengths.put(chr,chrs.get(chr).length);
			if(buildConsensus)
			{
				readSeqs.put(chr,new short[chrs.get(chr).length*4]);
			}
			if(verbose)System.err.print("Loaded "+chr+" consisting of "+chrs.get(chr).length+" nts. Mem total in MB: "+getGigsUsed()+"             \r");
		 }

		 fs.close();
		 System.gc();
		 System.err.println();
	}

	
	private void sortReads()
	{
		for(String str: reads.keySet())
		{
			Collections.sort(reads.get(str),this);
		}
		if(verbose)System.err.println("Finished sorting reads!");
	}
	
	private void processReads()
	{
		int totalReadsProcessed = 0;
		LinkedList<readWorker> workers = new LinkedList<readWorker>();
		LinkedList<String> chrStrs = new LinkedList<String>();
		int totalReads = 0;
		for(String str: reads.keySet())
		{
			chrStrs.add(str);
			totalReads+=reads.get(str).size();
		}
		int finished = 0;		                                
		if(!verbose) System.err.println("|0%                     25%                      50%                      75%                  100%|");
		while(finished < reads.size())
		{
			while(workers.size() < cores &&chrStrs.size() > 0)
			{
				String chrStr = chrStrs.removeFirst();
				workers.add(new readWorker(this,chrStr));
				if(verbose)System.err.println("Started processing "+chrStr+" with ~ "+reads.get(chrStr).size()/1000000+"M read regions. Total workers remaining: "+(reads.size()-finished));
				
			}
			LinkedList<readWorker> tempWorkers = new LinkedList<readWorker>();
			for(readWorker rw: workers)
			{
				if(!rw.isAlive())
				{
					finished++;
					totalReadsProcessed+= rw.count;
					double percentDone = ((100*((long)totalReadsProcessed))/totalReads);
					if(verbose)
						System.err.println("Finished processing "+rw.chrStr+" total % thus far: "+percentDone + " remaining workers: "+(reads.size()-finished)+" total memory: "+getGigsUsed()+" MB");
					else
					{
						String bar = "";
						for(int i = 0; i < percentDone; i++)
							bar+="=";
						System.err.print(bar+"                                               \r");
					}
				}
				else
					tempWorkers.add(rw);
			}
			workers = tempWorkers;
			tempWorkers = null;
			try {
				Thread.sleep(1000);  //wait for a sec
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		if(!verbose) System.err.println();
		//if(verbose)System.err.println("Finished building exons and defining junctions!");
	}
	
	
	private void outputExons() throws FileNotFoundException
	{
		if(verbose)System.err.println("Printing exons to file "+prefix+".exon.gff");
        PrintWriter pw = new PrintWriter(prefix+".exon.gff");
        pw.println("track name=\"MAPS Exons run on "+prefix+"\" description=\"MAPS predictions for "+prefix+"\" useScore=1");
        for(String chrStr: exons.keySet())
        {
        	String chr = chrStr.substring(0,chrStr.length()-1);
        	int count = 0;
        	for(Exon e: exons.get(chrStr))
        	{
        		count++;
                pw.println(chr+"\tMAPS\texon\t"+Math.max(1,e.start)+"\t"+(e.end-1)+"\t"+(e.getScore())+"\t+\t.\tE_"+chrStr+"_"+count);
        	}
        }
        pw.close();
	}
	
    private void outputJunctions() throws FileNotFoundException
    {
    	if(verbose)System.err.println("Printing junctions: "+junctions.size()+" to "+prefix+".junc.gff");
            PrintWriter pw = new PrintWriter(prefix+".junc.gff");
            pw.println("track name=\"MAPS Junctions run on "+prefix+"\" description=\"MAPS predictions for "+prefix+"\" useScore=1");
            int counter = 0;
            for(String start: junctions.keySet())
            {
                    HashMap<Integer,Integer> to = junctions.get(start);
                    
                    int best = 0;
                    int where = 0;
                    for(int end: to.keySet())
                    {
                            if(to.get(end) > best)
                            {
                                    best = to.get(end);
                                    where = end;
                            }
                    }
                    String chr = start.split("[+-]")[0];
                    int startPos = Math.max(1,Integer.parseInt(start.split("[+-]")[1]));
                    if(startPos > where)
                    {
                            int temp = where;
                            where = startPos;
                            startPos = temp;
                    }
                    int score = (int) Math.min(1000,10*best);

                    if(start.contains("+"))
                            pw.println(chr+"\tMAPS\tjunction\t"+startPos+"\t"+where+"\t"+score+"\t+\t.\tJ_"+chr+"+_"+counter);
                    else
                            pw.println(chr+"\tMAPS\tjunction\t"+startPos+"\t"+(where-1)+"\t"+score+"\t-\t.\tJ_"+chr+"-_"+counter);
                    
                    counter++;
            }
            pw.close();
    }	
	
	private void outputTranscripts(LinkedList<Transcript> transcripts, String chrStr, boolean plusStrand, PrintWriter fasta, PrintWriter gtf, PrintWriter full_gtf, PrintWriter nucFasta)
	{
		int tCount = 0;
		if(transcripts == null || transcripts.size()==0)
		{
			if(verbose) System.err.println(chrStr+" had no transcripts assembled. Skipping output.");
			return;
			
		}
		String chr = chrStr.substring(0,chrStr.length()-1);
		int[] transcriptSize = new int[10];
		int outputedFrags = 0;
		if(verbose)System.err.println("Outputing transcripts for: "+chrStr+" total memory: "+getGigsUsed());
		
		for(Transcript transcript: transcripts)
		{	
			boolean plus = plusStrand;
			//transcript.trimExtendedReads(ext); not required. We only use ext for filling gaps between reads... not extending them
			int transcriptStart = transcript.getStart();
			int transcriptEnd = transcript.getEnd();
			//end of exon search loop
			//lets output the transcript if it is of decent quality
			int exon_count = transcript.getExonCount();
			String transcriptSeq = "";
			String origSeq = "";
			for(Exon e: transcript.getExons())
			{
				transcriptSeq+=getSequence(chr,e.start-1,e.end-1,plus);
				origSeq+=getSeqWithoutConsensus(chr, e.start-1, e.end-1, plus);
			}	
			if(transcriptSeq.length() > 3)
			{
				tCount++;
				//lets output both positive and (if unstranded) the negative of the transcript
				String fullGTFp = "";
				String fullGTFm = "";
				String orfGTFp = "";
				String orfGTFm = "";
				String nucFASTAp = "";
				String pepFASTAp = "";
				String nucFASTAm = "";
				String pepFASTAm = "";
				int maxORFLength = 0;
				boolean plusLongest = true;
				for(int flip = 0; flip < 2; flip++)
				{
					if(unstranded && !plus)
					{
						//lets flip the transcript exons
						transcript.flipExons();
					}
					transcriptSize[Math.min(transcriptSize.length-1,exon_count-1)]++;

					//ok now lets output the entire transcript to GTF
					if(plus)
					{
						String tid ="MAPS."+chr+"."+tCount+"_+"+chr+":"+transcript.getStart()+"-"+transcript.getEnd();
						fullGTFp+=chr+"\tMAPS\ttranscript\t"+transcript.getStart()+"\t"+transcript.getEnd()+"\t"+transcript.getScore()+"\t+\t.\ttranscript_id \""+tid+"\"; gene_id \""+tid+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();
						for(Exon e: transcript.getExons())
						{
							fullGTFp+=chr+"\tMAPS\texon\t"+Math.max(1,e.start)+"\t"+(e.end-1)+"\t"+e.getScore()+"\t+\t.\ttranscript_id \""+tid+"\"; gene_id \""+tid+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\"; depth \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments,(avgFragLength*e.count)/(e.end-e.start))+System.lineSeparator();
						}
					}
					else
					{
						String tid ="MAPS."+chr+"."+tCount+"_-"+chr+":"+transcript.getStart()+"-"+transcript.getEnd();
						fullGTFm+=chr+"\tMAPS\ttranscript\t"+transcript.getStart()+"\t"+transcript.getEnd()+"\t"+transcript.getScore()+"\t-\t.\ttranscript_id \""+tid+"\"; gene_id \""+tid+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();
						for(Exon e: transcript.getExons())
						{
							fullGTFm+=chr+"\tMAPS\texon\t"+Math.max(1,e.start)+"\t"+(e.end-1)+"\t"+e.getScore()+"\t-\t.\ttranscript_id \""+tid+"\"; gene_id \""+tid+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\"; depth \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments,(avgFragLength*e.count)/(e.end-e.start))+System.lineSeparator();
						}
					}
					if(ap.get("seq") != null)  //we provided a genomic sequence so we have ORF-level assembly
					{
						for(int frame = 0; frame < 3; frame++)  //for each frame
						{	
							int preset = 0;
							int afterset = 0;
							String untranslated = transcriptSeq.substring(frame);
							int tail = untranslated.length()%3;
							if(tail > 0)
								untranslated = untranslated.substring(0, untranslated.length()-tail);
							String translated = seqToAA(untranslated);
							String origTrans = seqToAA(origSeq.substring(frame));							
							int actualFrame = (transcript.getStart()+frame)%3;
							if(!plus) actualFrame = (transcript.getEnd()-frame)%3;
							boolean appending = true;
							boolean oappending = true;
							while(extendTranscripts)
							{
								//lets go upstream until we hit an X
								String peek = "";
								String opeek = "";
								if(plus && transcriptStart+frame-preset-4 >= 0 && preset < 30000)
								{
									String otrip = getSeqWithoutConsensus(chr,transcriptStart+frame-preset-4,transcriptStart+frame-preset-1,plus);
									String trip = getSequence(chr,transcriptStart+frame-preset-4,transcriptStart+frame-preset-1,plus);
									peek = seqToAA(trip);
									opeek = seqToAA(otrip);
									if(appending && peek.compareTo("X")!= 0 && peek.compareTo("!")!= 0)
									{
										translated=peek+translated;
										untranslated=trip+untranslated;
										if(!oappending)
										{
											if(!origTrans.contains("@"))
												origTrans="@X"+origTrans;
											else
												origTrans="@"+origTrans;
										}
									}
									else 
									{
										appending =false;
										if(!oappending) break;
									}
									if(oappending && opeek.compareTo("X")!= 0 && opeek.compareTo("!")!= 0)
									{
										origTrans=opeek+origTrans;
										if(!appending)
										{
											if(!translated.contains("@"))
												translated="@X"+translated;
											else
												translated="@"+translated;
										}
									}
									else 
									{
										oappending =false;
										if(!appending) break;
									}
									preset+=3;
								}
								else if(!plus && transcriptEnd-frame+preset+2 < chrLengths.get(chr) && preset < 30000)
								{
									String trip =getSequence(chr,transcriptEnd-frame+preset-1,transcriptEnd-frame+preset+2,plus);
									String otrip =getSeqWithoutConsensus(chr,transcriptEnd-frame+preset-1,transcriptEnd-frame+preset+2,plus);
									peek = seqToAA(trip);
									opeek=seqToAA(otrip);
									if(appending && peek.compareTo("X")!= 0 && peek.compareTo("!")!= 0)
									{
										translated=peek+translated;
										untranslated = trip+untranslated;
										if(!oappending)
										{
											if(!origTrans.contains("@"))
												origTrans="@X"+origTrans;
											else
												origTrans="@"+origTrans;
										}
									}
									else 
									{
										appending = false;
										if(!oappending) break;
									}
									if(oappending && opeek.compareTo("X")!= 0 && opeek.compareTo("!")!= 0)
									{
										origTrans=opeek+origTrans;
										if(!appending)
										{
											if(!translated.contains("@"))
												translated="@X"+translated;
											else
												translated="@"+translated;
										}
									}
									else 
									{
										oappending = false;
										if(!appending) break;
									}
									preset+=3;
								}
								else
								{
									break;
								}
							}  //done with prepend
							appending = true;
							oappending =true;
							while(extendTranscripts)
							{
								
								String peek = "";
								String opeek = "";
								if(plus && transcriptEnd+frame+afterset+1 < chrLengths.get(chr) && afterset < 30000)
								{

									String trip = getSequence(chr,transcriptEnd+frame+afterset-2,transcriptEnd+frame+afterset+1,plus);
									String otrip = getSeqWithoutConsensus(chr,transcriptEnd+frame+afterset-2,transcriptEnd+frame+afterset+1,plus);
									peek = seqToAA(trip);
									opeek = seqToAA(otrip);
									if(appending && peek.compareTo("X")!= 0 && peek.compareTo("!")!= 0)
									{
										translated+=peek;
										untranslated+=trip;
										if(!oappending)
										{
											if(!origTrans.contains("@"))
												origTrans+="X@";
											else
												origTrans+="@";
										}
									}
									else
									{
										if(peek.compareTo("X")==0)
											origTrans+=opeek;
										appending = false;
										if(!oappending) break;
									}
									if(oappending && opeek.compareTo("X")!= 0 && opeek.compareTo("!")!= 0)
									{
										origTrans+=opeek;
										if(!appending)
										{
											if(!translated.contains("@"))
												translated+="X@";
											else
												translated+="@";
										}
									}
									else
									{
										oappending = false;
										if(!appending) break;
									}
									afterset+=3;
								}
								else if(!plus && transcriptStart-frame-afterset-2 >= 0 && afterset <30000)
								{
									String trip =getSequence(chr,transcriptStart-frame-afterset-2,transcriptStart-frame-afterset+1,plus);
									String otrip =getSeqWithoutConsensus(chr,transcriptStart-frame-afterset-2,transcriptStart-frame-afterset+1,plus);
									peek = seqToAA(trip);
									opeek = seqToAA(otrip);
									if(appending && peek.compareTo("X")!= 0 && peek.compareTo("!")!= 0)
									{
										translated+=peek;
										untranslated+=trip;
										if(!oappending)
										{
											if(!origTrans.contains("@"))
												origTrans+="X@";
											else
												origTrans+="@";
										}
									}
									else 
									{
										if(peek.compareTo("X")==0)
											origTrans+=opeek;
										appending = false;
										if(!oappending) break;
									}
									if(oappending && opeek.compareTo("X")!= 0 && opeek.compareTo("!")!= 0)
									{
										origTrans+=opeek;
										if(!appending)
										{
											if(!translated.contains("@"))
												translated+="X@";
											else
												translated+="@";
										}
									}
									else 
									{
										oappending = false;
										if(!appending) break;
									}
									afterset+=3;
								}
								else
								{
									break;
								}
							}//finished extending the sequence and the unmutated sequence
							String[] split = translated.split("X"); // break it up by stops to get the ORFs
							if(origTrans.length() == translated.length())
								origTrans+="X";
							String origTranslated = new String(origTrans);
							int fragCount = 0;
							int currentPos=0;
							if(plus)
								nucFASTAp+=">";
							else
								nucFASTAm+=">";
							for(String frag: split)
							{
								String tags = "";
								String nonCon= "";
								nonCon = origTrans.substring(0, Math.min(origTrans.length(),frag.length()+1));
								origTrans=origTrans.substring(Math.min(origTrans.length(),frag.length()+1));
								//lets add the mutation tags
								boolean shortened = false,elongated = false,nonSynonymous = false,removedStop = false;
								for(int c = 0; c < Math.min(nonCon.length(),frag.length()); c++)
								{
									if(frag.charAt(c) != nonCon.charAt(c))
									{
										//we had a difference
										if(frag.charAt(c)=='@')
											shortened = true;
										if(nonCon.charAt(c)=='@')
											elongated=true;
										if(nonCon.charAt(c)=='X')
											removedStop =true;
										else nonSynonymous=true;
									}
								}
								
								if(nonCon.length() > 0 && nonCon.charAt(nonCon.length()-1)!='X')
									tags = "[X]";
								if(shortened) tags+="[S]";
								if(elongated) tags+="[L]"; 
								if(nonSynonymous) tags+="[M]"; 
								if(removedStop) tags+="[G]";
								frag = frag.replaceAll("@", "");
								if(tags.length() > 0) tags= "_"+tags;
								fragCount++;
								if(frag.length() <= fmax && frag.length() >= fmin)
								{
									String fragNuc = untranslated.substring(currentPos*3+frame,Math.min(currentPos*3+frame+frag.length()*3,untranslated.length()));
									//scan the first couple of codons for alternate starts
									//'ATG','CTG','GTG','TTG','AAG','ACG','AGG','ATA','ATC','ATT'
									for(int i = 0; i < Math.min(fragNuc.length()/10,12);i+=3)
									{
										int score = 0;
										if(fragNuc.charAt(i)== 'A') score++;
										if(fragNuc.charAt(i+1)== 'T') score++;
										if(fragNuc.charAt(i+2)== 'G') score++;
										if(score == 2)
											tags+="_"+fragNuc.charAt(i)+fragNuc.charAt(i+1)+fragNuc.charAt(i+2);
										else if(score == 3)
											tags+="_ATG";
									}
									if(plus)
									{
										int left = Math.max(1,transcript.getRealPos(currentPos*3+frame,preset,plus));
										currentPos+=frag.length()+1;
										int right = Math.max(2,transcript.getRealPos((currentPos-1)*3+frame,preset,plus))-1;
										if(left > right)
										{
											int temp = left;
											left = right;
											right = temp;
										}
										if(right -left >maxORFLength)
										{
											maxORFLength = right-left;
											plusLongest=true;
										}
										String id ="MAPS."+chr+"."+tCount+"_+"+chr+":"+left+"-"+right+"_F:"+actualFrame+"_P:"+fragCount+tags;
										pepFASTAp+=">"+id+System.lineSeparator();;
										pepFASTAp+=frag+System.lineSeparator();
										nucFASTAp+=id+" ";
										outputedFrags++;
										//gtf output
										LinkedList<Exon> subset = transcript.getSubset(left,right,preset,afterset,plus);
										if(subset.size() == 0)
										{
											subset.add(transcript.get(0));
										}
										double avg = 0;
										for(Exon e: subset)
											avg+=e.getScore();
											
										if(subset.size() > 0)
											avg/=subset.size();
										else if(right < transcript.getStart())
											avg = transcript.getExons().getFirst().getScore();
										else
											avg = transcript.getExons().getLast().getScore();
										orfGTFp+=chr+"\tMAPS\ttranscript\t"+left+"\t"+right+"\t"+avg+"\t+\t.\ttranscript_id \""+id+"\"; gene_id \""+id+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();
										for(Exon e: subset)
										{
											int start = Math.max(1,e.start);
											int end = e.end-1;
											if(start > end)
											{
												int temp = start;
												start = end;
												end = temp;
											}
											orfGTFp+=chr+"\tMAPS\texon\t"+start+"\t"+end+"\t"+e.getScore()+"\t+\t.\ttranscript_id \""+id+"\"; gene_id \""+id+"\"; depth \""+String.format("%3.3f\";",(avgFragLength*e.count)/((e.end-e.start)))+String.format(" TPM \"%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();									
										}
		
									}
									else
									{
										int right = Math.max(1,transcript.getRealPos(currentPos*3+frame,preset,plus));
										currentPos+=frag.length()+1;
										int left = Math.max(1,transcript.getRealPos((currentPos-1)*3+frame,preset,plus));
										if(right -left >maxORFLength)
										{
											maxORFLength = right-left;
											plusLongest=false;
										}
										String id ="MAPS."+chr+"."+tCount+"_-"+chr+":"+left+"-"+(right-1)+"_F:"+actualFrame+"_P:"+fragCount+tags;
										pepFASTAm+=">"+id+System.lineSeparator();
										pepFASTAm+=frag+System.lineSeparator();
										nucFASTAm+=id+" ";

										outputedFrags++;
										//if(id.compareTo("MAPS.chr1.417_-chr1:29320-29529_F:1_P:11")==0)
										//	System.err.println();
										LinkedList<Exon> subset = transcript.getSubset(left,right,preset,afterset,plus);
										if(subset.size() == 0)
										{
											subset.add(transcript.get(0));
										}
										double avg = 0;
										for(Exon e: subset)
											avg+=e.getScore();
										if(subset.size() > 0)
											avg/=subset.size();
										else if(right < transcript.getStart())
											avg = transcript.getExons().getFirst().getScore();
										else
											avg = transcript.getExons().getLast().getScore();
										avg/=subset.size();
										orfGTFm+=chr+"\tMAPS\ttranscript\t"+left+"\t"+right+"\t"+avg+"\t-\t.\ttranscript_id \""+id+"\"; gene_id \""+id+"\"; TPM \""+String.format("%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();
										for(Exon e: subset)
										{
											int start = Math.max(1,e.start);
											int end = e.end-1;
											if(start > end)
											{
												int temp = start;
												start = end;
												end = temp;
											}
											orfGTFm+=chr+"\tMAPS\texon\t"+start+"\t"+end+"\t"+e.getScore()+"\t-\t.\ttranscript_id \""+id+"\"; gene_id \""+id+"\"; depth \""+String.format("%3.3f\";",(avgFragLength*e.count)/((e.end-e.start)))+String.format(" TPM \"%3.3f\"; fitQual \"%3.2f\"; FPKM \"%3.2f\";",transcript.getTPM(),transcript.getError(),(transcript.getFPKM()*1000000000.0)/totalFragments)+System.lineSeparator();								
										}
									}
								}
								else
									currentPos+=frag.length()+1;
								
							}
							if(plus)
							{
								nucFASTAp+=System.lineSeparator();
								nucFASTAp+=untranslated+System.lineSeparator();
							}
							else
							{
								nucFASTAm+=System.lineSeparator();
								nucFASTAm+=untranslated+System.lineSeparator();
							}
						}//frame
					}//seq check
					if(unstranded)
					{
						plus = !plus;
						transcriptSeq=revComplement(transcriptSeq);
					}
					else
						break;
				}//flip
				//ok now lets check which strand won out
				if(unstranded)
				{
					//lets pick the one with the longest ORF to print, otherwise we print it as is
					if(plusLongest)
					{
						if(fasta != null)
							fasta.print(pepFASTAp);
						if(nucFasta != null)
							nucFasta.print(nucFASTAp);
						if(full_gtf != null)
							full_gtf.print(fullGTFp);
						if(gtf != null)
							gtf.print(orfGTFp);
					}
					else
					{
						if(fasta != null)
							fasta.print(pepFASTAm);
						if(nucFasta != null)
							nucFasta.print(nucFASTAm);
						if(full_gtf != null)
							full_gtf.print(fullGTFm);
						if(gtf != null)
							gtf.print(orfGTFm);
					}
				}
				else
				{
					if(plusStrand)
					{
						if(fasta != null)
							fasta.print(pepFASTAp);
						if(nucFasta != null)
							nucFasta.print(nucFASTAp);
						if(full_gtf != null)
							full_gtf.print(fullGTFp);
						if(gtf != null)
							gtf.print(orfGTFp);
					}
					else
					{
						if(fasta != null)
							fasta.print(pepFASTAm);
						if(nucFasta != null)
							nucFasta.print(nucFASTAm);
						if(full_gtf != null)
							full_gtf.print(fullGTFm);
						if(gtf != null)
							gtf.print(orfGTFm);
					}
				}
				
			}//transcript length check
		}//transcript from transcriptList loop

		if(verbose)
		{
			if(ap.get("seq")!=null)
				System.err.println(String.format("There was a total of %d transcripts resulting in %d potential ORFs of size %d-%d aa.",tCount,outputedFrags,fmin,fmax));
			else
				System.err.println(String.format("There was a total of %d transcripts.",tCount));
		}
		if(totalORFs.get(chr) == null) totalORFs.put(chr, outputedFrags);
		else totalORFs.put(chr, totalORFs.get(chr)+outputedFrags);
	}

	private String getSequence(String chr, int start,
			int stop) {
		try{
			byte[] chars = new byte[(int) (stop-start)];
			if(chrs != null && chrs.get(chr) != null)
			{
				byte[] bytes = chrs.get(chr);
				try{
					for(int i = start; i < stop; i++)
					{
						chars[i-start]=bytes[i];
					}
				}catch(Exception e)
				{ return "!";}
				return new String(chars);
			}
			else
				return "!";
		}catch(Exception e)
		{
			return "!";
		}
	}
	
	private String getSequence(String chr, int start, int stop, boolean plus)
	{
		if(plus) return getSequence(chr,start,stop);
		else
			return revComplement(getSequence(chr,start,stop));
	}
	

	private String seqToAA(String seq)
	{		
		char[] arr = seq.toCharArray();
		char[] res = new char[arr.length/3];
		for(int i = 0; i <= seq.length()-3; i+=3 )
		{
			String trip = ""+arr[i]+arr[i+1]+arr[i+2];
			char aa = '!';
			if(nucToPep.get(trip) != null)
				aa = nucToPep.get(trip).charAt(0);
			res[i/3]=aa;
		}
		return new String(res);
	}
	
	private String getSeqWithoutConsensus(String chr, int start, int stop, boolean plus)
	{
		char[] nuc = getSequence(chr,start,stop).toCharArray();
		//now lets go through and find SNPs
		if(snps.get(chr) != null)
		{
			for(SNP snp: snps.get(chr))
			{
				if(snp.pos >= start && snp.pos < stop && nuc.length > snp.pos-start)
				{
					nuc[snp.pos-start]=snp.orig;
				}
			}
		}
		String unmutated = new String(nuc);
		if(!plus) unmutated = revComplement(unmutated);
		return unmutated;
	}
	
	private HashMap<String, String> buildNucToPep() {
		HashMap<String,String> map = new HashMap<String,String>();
		map.put("TTT","F");map.put("TCT","S");map.put("TAT","Y");map.put("TGT","C");
		map.put("TTC","F");map.put("TCC","S");map.put("TAC","Y");map.put("TGC","C");
		map.put("TTA","L");map.put("TCA","S");map.put("TAA","X");map.put("TGA","X");
		map.put("TTG","L");map.put("TCG","S");map.put("TAG","X");map.put("TGG","W");
		map.put("CTT","L");map.put("CCT","P");map.put("CAT","H");map.put("CGT","R");
		map.put("CTC","L");map.put("CCC","P");map.put("CAC","H");map.put("CGC","R");
		map.put("CTA","L");map.put("CCA","P");map.put("CAA","Q");map.put("CGA","R");
		map.put("CTG","L");map.put("CCG","P");map.put("CAG","Q");map.put("CGG","R");
		map.put("ATT","I");map.put("ACT","T");map.put("AAT","N");map.put("AGT","S");
		map.put("ATC","I");map.put("ACC","T");map.put("AAC","N");map.put("AGC","S");
		map.put("ATA","I");map.put("ACA","T");map.put("AAA","K");map.put("AGA","R");
		map.put("ATG","M");map.put("ACG","T");map.put("AAG","K");map.put("AGG","R");
		map.put("GTT","V");map.put("GCT","A");map.put("GAT","D");map.put("GGT","G");
		map.put("GTC","V");map.put("GCC","A");map.put("GAC","D");map.put("GGC","G");
		map.put("GTA","V");map.put("GCA","A");map.put("GAA","E");map.put("GGA","G");
		map.put("GTG","V");map.put("GCG","A");map.put("GAG","E");map.put("GGG","G");
		map.put("NNN","!");
		HashMap<String,String> newMap = new HashMap<String,String>();
		for(String key: map.keySet())
		{
			String copy = new String(key);
			String TtoU=copy.replace('T', 'U');
			newMap.put(TtoU,map.get(key));
			newMap.put(key,map.get(key));
			newMap.put(key.toLowerCase(),map.get(key).toLowerCase());
			newMap.put(TtoU.toLowerCase(), map.get(key).toLowerCase());
		}
		return newMap;	
	}
	

	/**
	 * Return END,END/EXON_END,NEXT_EXON_START
	 * @param cigar
	 * @return
	 */
	private int[] getReadRegions(String cigar, int start)
	{
		int[] result = new int[2];
		int s = start;
		int e = s;
		Matcher m = p.matcher(cigar);
		while(m.find())
		{
			String g = m.group(1);
			if(g.endsWith("D"))  //delete the genome (represents a gap in the read alignment) 
				e-=Integer.parseInt(g.substring(0,g.length()-1));
			//else if(g.endsWith("I")) //insert  (represent a loop in the read)
			//	e+=Integer.parseInt(g.substring(0,g.length()-1));
			else if(g.endsWith("N"))  //junc
			{
				result[result.length-1]=e; //end
				result[result.length-2]=s; //start
				int[] newResult= new int[result.length+2];
				for(int i = 0; i < result.length; i++)
					newResult[i]=result[i];
				result=newResult;
				e+=Integer.parseInt(g.substring(0,g.length()-1));
				s=e;
			}
			else if(g.endsWith("S")||g.endsWith("H"))  //assume that these only occur at the beginning or end of a read...
			{
				//do nothing (according to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/figure/F1/, SAM positions ignore these)
			}
			else if(g.endsWith("M"))  //match of some kind
				e+=Integer.parseInt(g.substring(0,g.length()-1));
		}
		result[result.length-1]=e;
		result[result.length-2]=s;
		return result;
	}
	
	private String getEditedReadSequence(String cigar, String seq)
	{
		String result = "";
		Matcher m = p.matcher(cigar);
		while(m.find())
		{
			String g = m.group(1);
			if(g.endsWith("D"))
			{
				int e = Integer.parseInt(g.substring(0,g.length()-1));
				for(int i = 0; i < e; i++)
					result+="N";
			}
			else if(g.endsWith("N"))
			{
				// do nothing
			}
			else if(g.endsWith("I"))
			{
				//insertion means we have extra sequence that is not annotated... need to trim it to ensure it is not read as a series of SNPs
				seq = seq.substring(Integer.parseInt(g.substring(0,g.length()-1)),seq.length());
			}
			else if(g.endsWith("S")||g.endsWith("H"))  //assume that these only occur at the beginning or end of a read...
			{
				seq = seq.substring(Integer.parseInt(g.substring(0,g.length()-1)),seq.length());
			}
			else if(g.endsWith("M"))
			{
				result+= seq.substring(0,Integer.parseInt(g.substring(0,g.length()-1)));
				seq = seq.substring(Integer.parseInt(g.substring(0,g.length()-1)),seq.length());
			}
		}
		return result;
	}
	
	
	private LinkedList<int[]> combineReadRegions(int[] read1, int[] read2)
	{
		ArrayList<Integer> a = new ArrayList<Integer>();
		ArrayList<Integer> b = new ArrayList<Integer>();
		for(int i: read1)
			a.add(i);
		for(int i: read2)
			b.add(i);
		//Collections.sort(a);  //we are guaranteed that this is the case by the cigar string construction...
		//Collections.sort(b);
		for(int i = 0; i < a.size(); i+=2)
		{
			int min = a.get(i);
			int max = a.get(i+1);
			
			for(int j = 0; j < b.size(); j+=2)
			{
				int min2 = b.get(j);
				int max2 = b.get(j+1);
				if(min2 <= max && max2 >= min)  //overlap
				{
					a.set(i, Math.min(min,min2));
					a.set(i+1, Math.max(max, max2));
					b.remove(j);
					b.remove(j);
					min = a.get(i);
					max = a.get(i+1);
					//merged overlapping region (decreasing b)
				}
			}
		}
		
		//The remaining regions should be returned as a separate list because they represent a separate piece of the transcript that may have "hidden" gaps
		LinkedList<int[]> list = new LinkedList<int[]>();
		int[] result = new int[a.size()];
		for(int i = 0; i < a.size(); i++)
		{
			result[i]=a.get(i);
		}
		list.add(result);
		if(b.size() > 0)
		{
			int[] result2 = new int[b.size()];
			for(int i = 0; i < b.size(); i++)
			{
				result2[i]=b.get(i);
			}
			list.add(result2);
		}
		return list;
	}

	public static void main(String[] args) {
			new MAPS(args);
	}
	
	private static String revComplement(String seq)
	{
		char[] arr = seq.toCharArray();
		char[] ret = new char[arr.length];
		for(int i = arr.length-1; i >= 0; i--)
 		{
			if(arr[i]=='A')ret[ret.length-i-1]='T';
			else if(arr[i]=='T')ret[ret.length-i-1]='A';
			else if(arr[i]=='C')ret[ret.length-i-1]='G';
			else if(arr[i]=='G')ret[ret.length-i-1]='C';
			else if(arr[i]=='N')ret[ret.length-i-1]='N';
		}
		return new String(ret);
	}

	@Override
	public int compare(int[] o1, int[] o2) {
		// SORT by minimum of the read 
		if(o1[0] > o2[0]) return 1;
		else if(o1[0] < o2[0]) return -1;
		else return 0;
	}
	
	public int getGigsUsed()
	{
		return (int) ((Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024));
	}

}
