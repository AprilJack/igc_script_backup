import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;


public class TranscriptAnnotater {
	/*
	 * -dataDir "/Volumes/viper.salk.edu/programs/TranscriptAnnotater/quick" -i "/Volumes/viper.salk.edu/Jiao/STARsep.gtf" -3frame "/Volumes/viper.salk.edu/Jiao/alignments/STARsep.fa" -only .gwas,.var,.bedGraph -filterannotated false -maxLength 10000 -minConservation 0
	 */

	HashMap<String,LinkedList<Transcript>> transcripts = new HashMap<String,LinkedList<Transcript>>();
	HashMap<String,LinkedList<Transcript>> tooLong = new HashMap<String,LinkedList<Transcript>>();
	HashMap<String,LinkedList<Transcript>> notConserved = new HashMap<String,LinkedList<Transcript>>();
	HashMap<String,LinkedList<Transcript>> tooannotated = new HashMap<String,LinkedList<Transcript>>();
	HashMap<String,Transcript> byId = new HashMap<String,Transcript>();
	
	LinkedList<String> allFiles = new LinkedList<String>();
	LinkedList<String> only = new LinkedList<String>();
	boolean stranded = false;
	boolean filter = true;
	int binSize = 1000000; 
	
	HashMap<String,HashMap<String,Annotation>> annotations = new HashMap<String,HashMap<String,Annotation>>();
	

	public TranscriptAnnotater(String[] args) throws IOException
	{
		ArgParser ag = new ArgParser("TranscriptAnnotater takes a gtf file of transcripts, and outputs the average expression, GWAS, Variant, and Conservation hits across exons.\n"
				+ "You can also first filter the transcripts by overlap with known anntoated transcripts by supplying a -reference GTF file. Gzipped (.gz) accepted.  Output tab delimited to stdOut.\n"
				+ "Assumes that your bedGraph files are normalized already. Conservation/Variant/GWAS annotation from UCSC accepted currently.");
		ag.registerArg("i", "input.gtf", "A gtf file containing unique transcripts with exons (required)");
		ag.registerArg("reference", "ref.gtf", "A gtf file describing the reference annotated transcripts and cds sites (required)");
		ag.registerArg("only", ".bedGraph,.gwas,.var", "list file substrings for datasets you want to use (bedGraph/gwas/uniprotvars). Datasets not containing any of these will be ommited.");
		ag.registerArg("conservation", "hg19100way.cons.gz", "Name of the conservation track with respect to dataDir. This assumes phastCons one-value per base.");
		ag.registerArg("dataDir", ".", "Directory containing other directories with bedGraph/GWAS/Conservation tracks. Searching is recursive.");
		ag.registerArg("maxLength", "10000", "transcript length must be less than this number of bp");
		ag.registerArg("filter", "false", "Remove transcripts with regions overlapping annotated regions");
		ag.registerArg("minConservation", "0", "transcripts are filtered by average conservation across exons.");
		ag.registerArg("stranded", "false", "bedGraph files are expected to have a positive and negative strand. Otherwise all are quantified.");	
		ag.parseArgs(args);
		if(args.length < 1)
		{
			ag.printUsage();
		}
		System.err.println(ag.getAllArgs().toString());
		String[] inc = ag.get("only").split(",");
		for(String s: inc)
		{
			only.add(s);
		}
		System.err.println("Will only search: "+ag.get("only"));
		stranded = ag.getAsBoolean("stranded");
		inputFile = ag.get("i");
		System.err.println("Loading "+inputFile);
		File input = new File(ag.get("i"));
		Scanner s =null;
		//read in the input GTF file containing ORFs or Transcripts
		if(input.getName().endsWith(".gz"))
			s = new Scanner(new GZIPInputStream(new FileInputStream(input)));
		else
			s = new Scanner(input);
		filter = ag.getAsBoolean("filter");
		if(filter) System.err.println("Will filter out transcripts overlapping with annotated expressed genes in "+ag.get("reference"));
		int totalTranscript = 0;
		long tic = System.currentTimeMillis();
		while(s.hasNextLine())
		{
			String line = s.nextLine();
			String[] split = line.split("\t");
			if(split.length > 8 && split[2].length() >0)
			{
				String chrStr = split[0];
				String ann = split[2];
				if(ann.compareTo("exon")==0)
				{
					int start = Integer.parseInt(split[3].trim());
					int end = Integer.parseInt(split[4].trim());
					int pos = start/binSize;
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));					
					if(byId.get(id) == null)
					{
						Transcript newT = new Transcript(id,chrStr,split[6].compareTo("+")==0);
						newT.addNote(split[1]); newT.addNote(split[8]);
						byId.put(id, newT);
						if(transcripts.get(chrStr+","+pos) == null) transcripts.put(chrStr+","+pos, new LinkedList<Transcript>());
						transcripts.get(chrStr+","+pos).add(byId.get(id));
						totalTranscript++;
					}
					byId.get(id).addRegion(start, end);


					if(totalTranscript %100000 == 0)System.err.print(".");
				}
			}
		}
		s.close();
		System.err.println();
		System.err.println("Loaded in "+byId.size()+" unique input transcripts");
		int totalCoding = 0;
		int totalAnnotations = 0;
		File refGTF = new File(ag.get("reference"));
		HashMap<String,Annotation> annId = new HashMap<String,Annotation>();
		if(refGTF.exists())
		{
			System.err.println("Loading "+ag.get("reference"));
			s =null;
			if(refGTF.getName().endsWith(".gz"))
				s = new Scanner(new GZIPInputStream(new FileInputStream(refGTF)));
			else
				s = new Scanner(refGTF);
			while(s.hasNextLine())
			{
				
				String line = s.nextLine();
				String[] split = line.split("\t");
				if(split.length > 8 && split[2].length() >0)
				{
					String chrStr = split[0];
					String ann = split[2];
					if(ann.compareTo("exon")==0)
					{
						int start = Integer.parseInt(split[3].trim());
						int end = Integer.parseInt(split[4].trim());
						int pos = start/binSize;
						String id = split[8].replaceAll("\"", "");
						id=id.replaceAll("\"", "");
						id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
						if(annId.get(id) == null)
						{
							if(annotations.get(chrStr+","+pos) == null)
								annotations.put(chrStr+","+pos, new HashMap<String,Annotation>());
							if(annotations.get(chrStr+","+pos).get(id) == null)
							{
								annotations.get(chrStr+","+pos).put(id, new Annotation(chrStr,id));
								totalAnnotations++;
							}
							
							annId.put(id,annotations.get(chrStr+","+pos).get(id));
						}
						annId.get(id).addRegion(start,end);
						
					}
					else if(ann.compareTo("CDS")==0)
					{
						//int start = Integer.parseInt(split[3].trim());
						//int end = Integer.parseInt(split[4].trim());
						String id = split[8].replaceAll("\"", "");
						id=id.replaceAll("\"", "");
						id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
						if(annId.get(id) != null && !annId.get(id).coding)
						{
							annId.get(id).coding=true;
							totalCoding++;
						}
					}
				}
				if(totalAnnotations % 10000 == 0)
					System.err.print(".");
			}
			s.close();
			System.err.println();
		}
		System.err.println("Loaded "+totalAnnotations+" annotated genes of which "+totalCoding+" were coding.");
		//lets filter the overlapping transcripts
		int transcriptsProcessed = 0;
		for(String chrStr: transcripts.keySet())
		{
			LinkedList<Transcript> remaining = new LinkedList<Transcript>();
			String[] split = chrStr.split(",");
			int pos = Integer.parseInt(split[1]);
			//now lets cross-reference against all nearby annotated genes
			LinkedList<Annotation> as = new LinkedList<Annotation>();
			if(annotations.get(chrStr) != null)
				as.addAll(annotations.get(chrStr).values());
			if(annotations.get(split[0]+","+(pos+1)) != null)
				as.addAll(annotations.get(split[0]+","+(pos+1)).values());
			for(Transcript t: transcripts.get(chrStr))
			{
				if(t.regions.size() ==0) continue;
				boolean hit = false;
				for(Annotation a: as)
				{
					if(a.coding && a.positions.size() > 0 && overlaps(a.getMin(), a.getMax(), t.getMin(), t.getMax()))
					{
						//lets search against the regions
						for(Region r: t.regions)
						{
							for(Region r2: a.positions)
							{
								if(r.overlaps(r2.start, r2.end)>0)
								{
									hit = true;
									break;
								}
							}
							if(hit) break;
						}
					}
				}
				if(hit)
				{
					//it is annotated
					t.annotated = true;
					t.similarity = 1;
					if(!filter) remaining.add(t);
				}
				else remaining.add(t);
				transcriptsProcessed++;
				if(transcriptsProcessed% 100000 == 0)
					System.err.print(".");
			}
			transcripts.put(chrStr,remaining);

		}
		
		int totalRemaining = 0;
		for(LinkedList<Transcript> list: transcripts.values())
		{
			totalRemaining+=list.size();
		}

		System.err.println();	
		System.err.println("Unannotated transcripts: "+totalRemaining+"/"+(totalTranscript)+" or about "+((100*totalRemaining)/(totalTranscript))+"%");
		System.err.println();
		int maxTranscriptLength = ag.getAsInt("maxLength");
		int totalOverall = 0;
		for(String chrStr: transcripts.keySet())
		{
			LinkedList<Transcript> remaining = new LinkedList<Transcript>();
			LinkedList<Transcript> longTranscripts = new LinkedList<Transcript>();
			int old = transcripts.get(chrStr).size();
			for(Transcript t: transcripts.get(chrStr))
			{
				totalOverall++;
				if(t.getTotalLength() <= maxTranscriptLength)
				{
					remaining.add(t);
				}
				else
					longTranscripts.add(t);
			}
			transcripts.put(chrStr,remaining);
			tooLong.put(chrStr,longTranscripts);
			//System.err.println(chrStr+" "+old+"->"+transcripts.get(chrStr).size());
		}
		System.err.println();
		if(filter)
			System.err.println("There were "+totalOverall+" unannotated transcripts that passed the size filter");
		else
			System.err.println("There were "+totalOverall+" transcripts that passed the size filter");

		File dataDir = new File(ag.get("dataDir"));
		//lets load the conservation next
		File f = new File(dataDir.getAbsolutePath()+File.separator+ag.get("conservation"));
		if(f.exists())
		{
			SuperScanner fs = new SuperScanner(f.getAbsolutePath());
			//this is a conservation track ... treat it similarly to bedGraph
			System.err.println("Loading conservation track: "+f.getName());
			//allFiles.add(f.getName());
			fs.getLine();  //header is ignored
			int count = 0;
			while(fs.hasMore())
			{
				count++;
				String line = fs.getLine();
				String[] split = line.split("\t");
				if(split.length > 5)
				{
					int cons = Integer.parseInt(split[5]);
					int start = Integer.parseInt(split[2]);
					int stop = Integer.parseInt(split[3]);
					String chr = split[1];
					int pos = start/binSize;
					if(transcripts.get(chr+","+pos) != null)
					{
						LinkedList<Transcript> ts = new LinkedList<Transcript>(); 
						if(transcripts.get(chr+","+pos) != null) ts.addAll(transcripts.get(chr+","+pos));
						if(transcripts.get(chr+","+(pos-1))!= null)
							ts.addAll(transcripts.get(chr+","+(pos-1)));
						for(Transcript t: ts)
						{
							if(t.overlaps(start, stop) >0)
								t.addConservation(f.getName(),cons,start,stop);
						}
					}
				}
				if(count %100000 ==0)
				{
					System.err.print(".");
				}
			}
			System.err.println();
		}
		//filter by conservation
		int conserved = 0;
		int minCons = ag.getAsInt("minConservation");
		for(String chr: transcripts.keySet())
		{
			LinkedList<Transcript> remaining = new LinkedList<Transcript>();
			LinkedList<Transcript> tooNew = new LinkedList<Transcript>();
			for(Transcript t: transcripts.get(chr))
			{
				if(t.conservation/t.getTotalLength() >= minCons)
				{
					conserved++;
					remaining.add(t);
				}
				else
					tooNew.add(t);
			}
			transcripts.put(chr,remaining);
			notConserved.put(chr,tooNew);
		}
		System.err.println("Conserved transcripts remaining: "+conserved+" (conservation score > "+minCons+")");
		
		if(dataDir.exists())
		{
			recursiveProcessing(dataDir);
		}
		
		System.err.println("Done loading all annotation data...");
		LinkedList<Transcript> sorted = new LinkedList<Transcript>();
		LinkedList<Transcript> bad = new LinkedList<Transcript>();
		for(String chr: transcripts.keySet())
		{
			for(Transcript t: transcripts.get(chr))
				sorted.add(t);
		}
		Collections.sort(sorted);
		for(String chr: tooLong.keySet())
		{
			for(Transcript t: tooLong.get(chr))
				bad.add(t);
		}
		for(String chr: tooannotated.keySet())
		{
			for(Transcript t: tooannotated.get(chr))
				bad.add(t);
		}
		for(String chr: notConserved.keySet())
		{
			for(Transcript t: notConserved.get(chr))
				bad.add(t);
		}
		sorted.addAll(bad);
		System.err.println("Total time: "+((System.currentTimeMillis()-tic)/60000)+" mins.");
		String header = "Transcripts "+(new Date()).toString()+"\tlength\tStrand\t#exons\toverlaps_ref?\tConservation\tChr\tStart\tEnd";
		for(String file: allFiles)
			header+="\t"+file;
		header+="\tNotes from "+ag.get("i");
		System.out.println(header);
		for(Transcript t: sorted)
		{
			StringBuilder sb = new StringBuilder(1000);
			sb.append(t.id);
			
			sb.append("\t");
			sb.append(t.getTotalLength());
			sb.append("\t");
			sb.append(t.regions.size());
			sb.append("\t");
			sb.append(t.plus ? "+" : "-");
			sb.append("\t");
			sb.append(t.chr);
			sb.append("\t");
			sb.append(t.regions.getFirst().start);
			sb.append("\t");
			sb.append(t.regions.getLast().end);
			sb.append("\t");
			sb.append((t.similarity > 0)? "TRUE" : "FALSE");
			sb.append("\t");
			sb.append(t.conservation/t.getTotalLength());
			
			for(String file: allFiles)
			{
				sb.append("\t");
				sb.append(t.toString(file));
			}
			sb.append("\t");
			for(String note: t.notes)
			{
				sb.append(note+" ");
			}
			System.out.println(sb.toString());
		}
	}
	

	private void recursiveProcessing(File dir)
	{
		for(File f: dir.listFiles())
		{
			if(f.isDirectory()) 
				recursiveProcessing(f);
			else
			{
				for(String keyword: only)
				{
					if(f.getName().contains(keyword) && f.getAbsolutePath().compareTo(inputFile)!= 0)
						annotateTranscripts(f);
				}
			}
		}
	}
	
	public void annotateTranscripts(File f)
	{
		long tic = System.currentTimeMillis();
		SuperScanner fs = null;
		try{

			fs = new SuperScanner(f.getAbsolutePath());
			if(f.getName().contains("bedGraph")){
				System.err.println("Loading bedGraph: "+f.getName());
				allFiles.add(f.getName());
				String signChar = "+";
				int count = 0;
				while(fs.hasMore())
				{
					count++;
					String line = fs.getLine();
	
					//we want to read this in as a bedGraph
					if(line.length()> 0)
					{
						if(line.startsWith("track"))
						{
							if(line.indexOf("+ strand")!= -1)
							{
								//we are processing the positive strand now
								signChar = "+";
							}
							else if(line.indexOf("- strand")!= -1)
							{
								//we are processing the positive strand now
								signChar = "-";
							}
						}
						else
						{
							String[] split = line.split("\t");
							if(split.length > 3)
							{
								float val = Float.parseFloat(split[3]);
								int start = Integer.parseInt(split[1]);
								int end = Integer.parseInt(split[2]);

								int pos = start/binSize;
								String chr = split[0];

								if(transcripts.get(chr+","+pos) != null)
								{
									LinkedList<Transcript> ts = new LinkedList<Transcript>(); 
									if(transcripts.get(chr+","+pos) != null) ts.addAll(transcripts.get(chr+","+pos));
									if(transcripts.get(chr+","+(pos-1))!= null)
										ts.addAll(transcripts.get(chr+","+(pos-1)));
									for(Transcript t: ts)
									{
										if(t.overlaps(start,end)>0 &&(!stranded||t.plus==(signChar.compareTo("+")==0)))
										{
											t.addValue(f.getName(),val,start,end);
										}
									}
								}
								
							}
						}
					}
					if(count % 1000000 == 0) System.err.print(".");
				}
				//now lets check each transcript and leave a 0 if it wasn't populated
				for(String chr: transcripts.keySet())
				{
					for(Transcript t: transcripts.get(chr))
					{
						if(t.values.get(f.getName())== null)
							t.values.put(f.getName(), new Value(0));
					}
				}
				System.err.println();
			}
			else if(f.getName().contains(".gwas"))
			{
				System.err.println("Loading GWAS variants file: "+f.getName());
				allFiles.add(f.getName());
				fs.getLine();  //header is ignored
				int count = 0;
				while(fs.hasMore())
				{
					count++;
					String line = fs.getLine();
					String[] split = line.split("\t");
					if(split.length > 17)
					{
						String annotation = split[4]+" "+split[5]+" "+split[10]+" "+split[17];
						int start = Integer.parseInt(split[2]);
						int stop = Integer.parseInt(split[3]);
						String chr = split[1];
						int pos = start/binSize;
						if(transcripts.get(chr+","+pos) != null)
						{
							LinkedList<Transcript> ts = new LinkedList<Transcript>(); if(transcripts.get(chr) != null) ts.addAll(transcripts.get(chr));
							if(transcripts.get(chr+","+(pos-1))!= null)
								ts.addAll(transcripts.get(chr+","+(pos-1)));
							for(Transcript t: ts)
							{
								if(t.overlaps(start,stop)>0 )
								{
									t.addValue(f.getName(),annotation,start,stop);
								}
							}
						}
					}
					if(count % 1000 ==0)
						System.err.print(".");
						
				}
				System.err.println();
			}
			else if(f.getName().contains(".vars"))
			{
				System.err.println("Loading UniprotVariants file: "+f.getName());
				allFiles.add(f.getName());
				fs.getLine();  //header is ignored
				int count = 0;
				while(fs.hasMore())
				{
					count++;
					String line = fs.getLine();
					String[] split = line.split("\t");
					if(split.length > 18)
					{
						String annotation = split[13]+" "+split[18];
						int start = Integer.parseInt(split[1]);
						int stop = Integer.parseInt(split[2]);
						String chr = split[0];
						int pos = start/binSize;
						if(transcripts.get(chr+","+pos) != null)
						{
							LinkedList<Transcript> ts = new LinkedList<Transcript>(); if(transcripts.get(chr) != null) ts.addAll(transcripts.get(chr));
							if(transcripts.get(chr+","+(pos-1))!= null)
								ts.addAll(transcripts.get(chr+","+(pos-1)));
							for(Transcript t: ts)							
							{
								if(t.overlaps(start,stop) > 0)
								{
									t.addValue(f.getName(),annotation,start,stop);
								}
							}
						}
					}
					if(count % 1000 == 0)
						System.err.print(".");
				}
				System.err.println();
			}
						
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally
		{
			fs.close();
		}
		System.err.println("Took "+((System.currentTimeMillis()-tic)/1000)+" seconds.");
	}
	


	private class Transcript implements Comparable<Transcript>
	{
		LinkedList<Region> regions = new LinkedList<Region>();
		String chr;
		String id;
		boolean plus;
		boolean annotated = false;
		double similarity = 0;
		private double conservation = 0;
		HashMap<String,Value> values = new HashMap<String,Value>();
		LinkedList<String> notes = new LinkedList<String>();
		
		
		public void addNote(String s)
		{
			notes.add(s);
		}
		
		public int getTotalLength()
		{
			int result = 0;
			for(Region r: regions)
			{
				result+=Math.abs(r.end-r.start);
			}
			return result;
		}
		
		public int getMin() {
			return Math.min(regions.getFirst().start,regions.getLast().start);
		}
		
		public int getMax() {
			return Math.max(regions.getFirst().end,regions.getLast().end);
		}

		public Transcript(String id,  String chrStr, boolean plus)
		{
			this.id=id; this.plus = plus; this.chr = chrStr;
		}
		
		public void addRegion(int start, int end)
		{
			int index =0; 
			for(Region r: regions)
			{
				if(r.start > start)
					break;
				index++;
			}
			regions.add(index,new Region(start,end));
		}
		
		public void addValue(String source,String val,int left, int right)
		{
			int overlap = overlaps(left,right);
			if(overlap > 0)
			{
				if(values.get(source) == null) values.put(source, new Value(val));
				else
					values.get(source).str+="; "+val;
			}
		}
		
		public void addConservation(String source,double val,int left, int right)
		{
			int overlap = overlaps(left,right);
			if(overlap > 0)
			{
				conservation+=val*overlap;
			}
		}
		
		public void addAnnotation(String source,String val)
		{

			if(values.get(source) == null) values.put(source, new Value(val));
			else
				values.get(source).str+="; "+val;
			annotated=true;
		}

		public void addValue(String source, double value,int left, int right)
		{
			int overlap =overlaps(left,right);
			if(values.get(source) == null) values.put(source, new Value(value*overlap));
			else
				values.get(source).val+=value*overlap;
		}
		
		public int overlaps(int left, int right)
		{
			if(regions.getFirst().start >= right || regions.getLast().end <= left)
				return 0;
			int total = 0;
			for(Region r: regions)
			{
				total+=Math.max(0,r.overlaps(left, right));
			}
			return total;
		}

		public Object toString(String file) {
			if(values.get(file) != null)
			{
				return values.get(file).toString(getTotalLength());
			}
			else return "";
		}


		@Override
		public int compareTo(Transcript t) {
			// TODO Auto-generated method stub
			double s1 = getScore();
			double s2 = t.getScore();
			if(s1 > s2) return -1;
			else if(s1 < s2) return 1;
			else return 0;
		}

		public double getScore() {
			// TODO Auto-generated method stub
			double score =conservation;
			if(annotated)
				score/=10;
			return score;
		}

		
	}
	
	private class Value
	{
		double val;
		String str;
		public Value(double val)
		{
			this.val = val;
		}
		
		public Value(String str)
		{
			this.str = str;
		}
		
		public String toString(int totalLength)
		{
			if(str != null)
				return str;
			else return String.format("%3.5f",Math.max(0.00001,Math.max(0.1,val)/totalLength));
		}
	}
	

	String inputFile = "";
	
	public static int overlap(int start, int end, int left, int right)
	{
		return Math.min(right,end)-Math.max(left,start);
	}
	
	public static boolean overlaps (int start, int end, int left, int right)
	{
		return overlap(start,end,left,right) > 0;
	}
	

	
	public static void main(String[] args) {
		try {
			new TranscriptAnnotater(args);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}

