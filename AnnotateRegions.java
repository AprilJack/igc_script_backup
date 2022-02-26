import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;


public class AnnotateRegions {

	LinkedList<String> regions = new LinkedList<String>();
	File gtf;
	File cons;
	HashMap<String,HashMap<String,LinkedList<GraphBlock>>> bedGraphs = new HashMap<String,HashMap<String,LinkedList<GraphBlock>>>();
	String regionFile = "";
	HashMap<String,LinkedList<Record>> records = new HashMap<String,LinkedList<Record>>();
	HashMap<String,LinkedList<Record>> records2 = new HashMap<String,LinkedList<Record>>();
	HashMap<String,LinkedList<ConBlock>> conservation = new HashMap<String,LinkedList<ConBlock>>();
	int maxCons = 0;
	HashMap<String,Integer> consCounts = new HashMap<String,Integer>(); //naive assumption that 99% of genome will be covered by all species
	
	private class Record implements Comparable<Record>
	{
		String annotation;
		String chr; 
		int start;
		int end;
		String plusString = "+";

		public Record(String a, int s, int e, String chr, boolean plus)
		{
			annotation = a; start = s; end = e;this.chr = chr;
			if(!plus) plusString = "-";
			//lets trim the annotation
		}

		@Override
		public int compareTo(Record o) {
			// TODO Auto-generated method stub
			return o.start-start;
		}
		
		public boolean overlaps(int left, int right, String chr,String plusString)
		{
			if(this.chr.compareTo(chr)!= 0 || this.plusString.compareTo(plusString)!=0) return false; //wrong chr
			if((left >=start && left <= end)||(right >=start &&right <= end)||(start >=left && start <= right)||(end >=left &&end <= right))
				return true;
			return false;
		}
	}
	
	private class GraphBlock
	{
		int start;
		int end;
		float val;
		public GraphBlock(int start, int end, float val)
		{
			this.start=start;this.end=end; this.val=val;
		}
	}
	
	private class ConBlock
	{
		int start;
		int end;
		int score;
		public ConBlock(int start, int end, int score)
		{
			this.start=start;this.end=end;this.score=score;
		}
	}
	
	
	public AnnotateRegions(String[] args)
	{
		ArgParser ag = new ArgParser("annotateRegions extracts regions for transcripts and outputs conservation and overlap with known genomic features\n"
				+ "Outputs annotated txt file to stdout");
		ag.registerArg("i", null, "list of regions to annotation using the \"chr:min-max\" format. If not specified then regions are prompted");
		ag.registerArg("gtf", null, "Genomic GTF file containing information about known gene annotations");
		ag.registerArg("cons", null, "table file with col 2 = chr, 3 = start, 4 = end, score = 6. First line header is removed.");
		ag.registerArg("bedGraphs", null, "bedGraph files separated by a comma. The average value across each region will be recorded in the output.");
		ag.parseArgs(args);
		if(args.length < 1)
		{
			ag.printUsage();
		}
		gtf = new File(ag.get("gtf"));
		cons = new File(ag.get("cons"));
		regionFile= ag.get("i");
		String[] bedGraphFiles = ag.get("bedGraphs").split(",");
		if(bedGraphFiles.length > 0)
		{
			for(int i = 0; i < bedGraphFiles.length; i++)
			{
				if((new File(bedGraphFiles[i])).exists())
				{
					FastScanner fs=  new FastScanner(bedGraphFiles[i]);
					HashMap<String,LinkedList<GraphBlock>> graph = new HashMap<String,LinkedList<GraphBlock>>();
					//int lastCoord = 0;
					String signChar = "+";
					while(fs.hasMore())
					{
						String line = fs.getLine();
						if(line.length()> 0)
						{
							if(line.indexOf("+ strand")> -1)
							{
								//we are processing the positive strand now
								signChar = "+";
							}
							else if(line.indexOf("- strand")> -1)
							{
								//we are processing the positive strand now
								signChar = "-";
							}
							else
							{
								String[] split = line.split("\t");
								if(graph.get(split[0]+signChar)==null)
								{
									graph.put(split[0]+signChar, new LinkedList<GraphBlock>()); 
								}
								float val = Float.parseFloat(split[3]);
								int start = Integer.parseInt(split[1]);
								int end = Integer.parseInt(split[2]);
								graph.get(split[0]+signChar).add(new GraphBlock(start,end,val));
							}
						}
						
					}
					bedGraphs.put(bedGraphFiles[i], graph);
					System.err.println("Loaded "+bedGraphFiles[i]);
				}
				else
					System.err.println("File "+bedGraphFiles[i]+" doesn't exist!");
			}
		}

		if(gtf.exists())
		{
			System.err.println("Loading "+gtf.getName());
			FastScanner fs = new FastScanner(gtf.getAbsolutePath());
			while(fs.hasMore())
			{
				String line = fs.getLine();
				String[] split = line.split("\t");
				if(split.length > 8 && split[2].length() >0)
				{
					String chrStr = split[0];
					String ann = split[2];
					int start = Integer.parseInt(split[3]);
					int end = Integer.parseInt(split[4]);
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					Record r = new Record(ann,start,end,chrStr,split[6].compareTo("+")==0);
					if(records.get(id)== null) records.put(id,new LinkedList<Record>());
					records.get(id).add(r);					
				}
			}
		}
		if(cons.exists())
		{
			System.err.println("Loading "+cons.getName());
			FastScanner fs = new FastScanner(cons.getAbsolutePath());
			fs.getLine();  //header is ignored
			while(fs.hasMore())
			{
				String line = fs.getLine();
				String[] split = line.split("\t");
				if(split.length > 5)
				{
					if(conservation.get(split[1])==null) conservation.put(split[1], new LinkedList<ConBlock>());
					//lets save it
					int cons = Integer.parseInt(split[5]);
					conservation.get(split[1]).add(new ConBlock(Integer.parseInt(split[2]),Integer.parseInt(split[3]),cons));
					if(cons > maxCons) maxCons = cons;
				}
			}
			fs.close();
		}
		String header = "Region";
		if(cons != null)
			header+="\t"+cons.getName();
		for( String graph: bedGraphs.keySet())
			header+="\t"+graph;
		if(gtf != null)
			header+="\t"+gtf.getName();
		System.out.println(header);
		if(regionFile != null)
		{
			System.err.println("Loading "+regionFile);
			FastScanner fs = new FastScanner(regionFile);
			//String currentTranscript = "";
			while(fs.hasMore())
			{
				String line = fs.getLine();
				String[] split = line.split("\t");
				if(split.length > 8 && split[2].length() >0)
				{
					String chrStr = split[0];
					String ann = split[2];
					int start = Integer.parseInt(split[3]);
					int end = Integer.parseInt(split[4]);
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					Record r = new Record(ann,start,end,chrStr,split[6].compareTo("+")==0);
					if(records2.get(id)== null) records2.put(id,new LinkedList<Record>());
					records2.get(id).add(r);					
				}
			}
		}

		Scanner s = new Scanner(System.in);
		while(true)
		{
			System.out.println("Please specify a transcript ID");
			String line = s.nextLine();
			try{
				if(line.length() > 0)
				{
					LinkedList<Record> transcript = records2.get(line);
					if(transcript != null)
					{
						outputRegion(line,transcript);
					}
					else
						System.out.println("Could not find "+line+" in the input file");
				}
				else
					break;
			}catch(Exception e){System.out.println("Bad format try again!");}
		}
		s.close();
		
	}
	
	private void outputRegion(String transcript,LinkedList<Record> list)
	{
		String output = transcript+"\t";
		String annotations = "";
		int totalLength = 0;
		double consScore = 0;
		HashMap<String,Double> bedGraphScores = new HashMap<String,Double>();
		for(Record r: list)
		{
			int min = r.start;
			int max = r.end;
			if(cons.exists())
			{

			//	int count = 0;
				try{
					LinkedList<ConBlock> blocks = conservation.get(r.chr);
					for(ConBlock b: blocks)
					{
						if(b.start<= max && b.end >= min)
						{
							//lets compute the overlap
							int length = Math.min(b.end,max)-Math.max(b.start,min);
							consScore+=length*(100.0*b.score)/maxCons;
							totalLength+=length;
						//	count+=length;
						}
					}
				}catch(Exception e)
				{
					System.err.println(e.getMessage());
				}
			}
			for(String graph: bedGraphs.keySet())
			{
				LinkedList<GraphBlock> blocks = bedGraphs.get(graph).get(r.chr+r.plusString);
				if(bedGraphScores.get(graph) == null) bedGraphScores.put(graph,0.0);
				for(GraphBlock gb : blocks)
				{
					if(gb.start<= max && gb.end >= min)
					{
						//lets compute the overlap
						int length = Math.min(gb.end,max)-Math.max(gb.start,min);
						bedGraphScores.put(graph,bedGraphScores.get(graph)+length*gb.val);
					}
				}				
			}
			if(gtf.exists())
			{
				for(String feature:records.keySet())
				{
					LinkedList<Record> recs = records.get(feature);
					for(Record r2: recs)
					{
						if(r2.overlaps(min, max,r.chr,r.plusString) && transcript.length() > 0 && r2.annotation.length() > 0)
						{
							annotations+=feature+" "+r2.annotation+";";
						//	break;
						}
					}
				}
			}
		}
		//now lets output the whole enchillada
		output+=String.format("%3.3f",consScore/totalLength);
		for(String graph: bedGraphScores.keySet())
		{
			output+=String.format("\t%3.3f",bedGraphScores.get(graph)/totalLength);
		}
		if(annotations.length() == 0)
			annotations="Unannotated";
		output+="\t"+annotations;	
		System.out.println(output);
	}
	
	public static void main(String[] args) {
		new AnnotateRegions(args);

	}

}
