import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;


public class annotateRegion {

	LinkedList<String> regions = new LinkedList<String>();
	File gtf;
	File cons;
	String regionFile = "";
	HashMap<String,LinkedList<Record>> records = new HashMap<String,LinkedList<Record>>();
	HashMap<String,double[]> conservation = new HashMap<String,double[]>();
	HashMap<String,Integer> consCounts = new HashMap<String,Integer>(); //naive assumption that 99% of genome will be covered by all species
	
	private class Record implements Comparable<Record>
	{
		String annotation;
		int start;
		int end;
		String id;
		
		public Record(String a, int s, int e, String i)
		{
			annotation = a; start = s; end = e; id = i;
		}

		@Override
		public int compareTo(Record o) {
			// TODO Auto-generated method stub
			return o.start-start;
		}
		
		public boolean overlaps(int left, int right)
		{
			if((left <= end && left >= start)||(right <= end && right >= start)||(start <= right && start >= left)||(end <= right && right >= start))
			{
				return true;
			}
			return false;
		}
	}
	
	
	
	public annotateRegion(String[] args)
	{
		ArgParser ag = new ArgParser("annotateRegionsExtracts regions for transcripts and outputs conservation and overlap with known genomic features\n"
				+ "Outputs annotated txt file to stdout");
		ag.registerArg("i", "selected.txt", "list of regions to annotation using the \"chr:min-max\" format. If not specified then regions are prompted");
		ag.registerArg("gtf", "genome.gtf", "Genomic GTF file containing information about known gene annotations");
		ag.registerArg("cons", "cons.wig", "fixed length wig file containing information about chr-specific conservation scores");
		ag.parseArgs(args);
		if(args.length < 1)
		{
			ag.printUsage();
		}
		gtf = new File(ag.get("gtf"));
		cons = new File(ag.get("cons"));
		regionFile= ag.get("i");
		if(gtf.exists())
		{
			System.err.println("Loading "+gtf.getName());
			FastScanner fs = new FastScanner(gtf.getAbsolutePath());
			while(fs.hasMore())
			{
				String line = fs.getLine();
				String[] split = line.split("\t");
				String chrStr = split[0]+split[6];
				String ann = split[2];
				int start = Integer.parseInt(split[3]);
				int end = Integer.parseInt(split[4]);
				String id = split[8].replaceAll("\"", "");
				Record r = new Record(ann,start,end,id);
				if(records.get(chrStr)== null) records.put(chrStr,new LinkedList<Record>());
				records.get(chrStr).add(r);
			}
		}
		if(cons.exists())
		{
			System.err.println("Loading "+cons.getName());
			FastScanner fs = new FastScanner(cons.getAbsolutePath());
			String chr = null;
			int start = 0;
			LinkedList<Double> values = new LinkedList<Double>();
			while(fs.hasMore())
			{
				String line = fs.getLine();
				if(line.startsWith("fixedStep"))
				{
					if(chr != null)
					{
						//lets save the chr
						double[] arr = new double[values.size()];
						int index=0;
						for(double d: values)
						{
							arr[index]=d;
							index++;
						}
						if(conservation.get(chr) == null)
						{
							conservation.put(chr, arr);
						}
						else
						{
							//lets add it
							double[] temp = conservation.get(chr);
							double[] combined = new double[Math.max(temp.length, arr.length)];
							for(int i = 0; i < combined.length; i++)
							{
								if(temp.length > i)
									combined[i]+=temp[i];
								if(arr.length > i)
									combined[i]+=arr[i];
							}
						}
						if(consCounts.get(chr) == null) consCounts.put(chr,0);
						consCounts.put(chr, consCounts.get(chr)+1);
						System.err.println("Loaded "+chr+" containing "+values.size());
					}
					chr = line.split("\\s+")[1].split("=")[1];
					start = Integer.parseInt(line.split("\\s+")[2].split("=")[1]);
					values = new LinkedList<Double>();
					for(int i = 0; i < start; i++)
						values.add(0.0);
					
				}
				else
					values.add(Double.parseDouble(line));
			}
			fs.close();
		}
		if(regionFile != null)
		{
			FastScanner fs = new FastScanner(regionFile);
			while(fs.hasMore())
			{
				String line = fs.getLine();
				//lets parse the region
				String[] split = line.split(":");
				String chr = split[0];
				String[] minmax = split[1].split("-");
				int min = Integer.parseInt(minmax[0]);
				int max = Integer.parseInt(minmax[1]);
				outputRegion(chr,min,max);
			}
			fs.close();
		}
		else
		{
			Scanner s = new Scanner(System.in);
			while(true)
			{
				String line = s.nextLine();
				if(line.length() > 0)
				{
					String[] split = line.split(":");
					String chr = split[0];
					String[] minmax = split[1].split("-");
					int min = Integer.parseInt(minmax[0]);
					int max = Integer.parseInt(minmax[1]);
					outputRegion(chr,min,max);
				}
				else
					break;
			}
			s.close();
		}
		
	}
	
	private void outputRegion(String chr, int min, int max)
	{
		String output = chr+":"+min+"-"+max;
		if(cons.exists())
		{
			double avg = 0;
			int count = 0;
			try{
				for(int i = min; i <= max; i++)
				{
					avg+=conservation.get(chr)[i];
					count++;
				}
			}catch(Exception e)
			{
				System.err.println(e.getMessage());
			}
			output+="\t"+(avg/(count*consCounts.get(chr)));
		}
		if(gtf.exists())
		{
			LinkedList<Record> recs = records.get(chr);
			boolean overlaps = false;
			for(Record r: recs)
			{
				if(r.overlaps(min, max))
				{
					output+="\tOverlaps";
					overlaps = true;
					break;
				}
			}
			if(!overlaps)
				output+="\tNo overlap";
		}
		System.out.println(output);
	}
	
	public static void main(String[] args) {
		new annotateRegion(args);

	}

}
