import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;


public class countTagsInRegions {

	int res = 100000;
	int dirs = 1;
	double threshold=0.1;
	double maxRepThresh = 1.0;
	HashMap<String,LinkedList<Integer>> expNames = new HashMap<String,LinkedList<Integer>>();
	
		//  chr           pos/res   counts
	HashMap<String,HashMap<String,int[]>> map = new HashMap<String,HashMap<String,int[]>>();
	
	public countTagsInRegions(String[] args) {
		ArgParser parser =new ArgParser("countTagsInRegions takes HIC PE Homer Tag directories and outputs a list of all interactions (regardless of distance) for processing with differential expression tools");
		parser.registerArg("res", "100000", "fixed window resolution to use for counting up interaction reads");
		parser.registerArg("threshold", "0.1", "only interactions that have (max-min)/(1+min+max) interactions are kept. This ensures variable interactions across conditions");
		parser.registerArg("exp", null, "Comma-separated exp names (e.g. WT,WT,MT,MT). Describes the experimental setup of your samples. Required if you want to ensure that only interactions with comparable replicate values are kept.");
		parser.registerArg("maxRepThresh", "1.0","Uses exp to figure out condition replicates and imposes a maximum (max-min)/(1+min+max) for valid interactions. I.e. keep interactions with relatively similar replicate values");
		parser.parseArgs(args);
		System.err.println(parser.getAllArgs());
		if(parser.getList().size() < 1)
		{
			parser.printUsage();
			System.exit(1);
		}
		dirs = parser.getList().size();
		res = parser.getAsInt("res");
		threshold = parser.getAsDouble("threshold");
		maxRepThresh = parser.getAsDouble("maxRepThresh");
		if(parser.get("exp")!= null)
		{
			String[] names = parser.get("exp").split(",");
			int index = 0;
			for(String name: names)
			{
				if(expNames.get(name) == null) expNames.put(name, new LinkedList<Integer>());
				expNames.get(name).add(index);
				index++;
			}
			if(index != dirs)
			{
				System.err.println("Number of exp names ("+index+") doesn't match the number of tag dirs specified: "+dirs);
				System.exit(1);
			}
		}
		long[] totalCounts = new long[dirs];
		LinkedList<String> dirNames= parser.getList();
		int i =0;
		for(String dirName : dirNames)
		{
			File dir = new File(dirName);
			for(File f: dir.listFiles())
			{
				if(f.getName().endsWith("tags.tsv"))
				{
					//lets load it in 
					try {
						int count = loadIn(f,i);
						System.err.println("Loaded "+dirName+": "+f.getName()+" with "+count);
						totalCounts[i]+=count;
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			i++;
			
		}
		//now lets output all of the interactions
		int count = 1;
		String header = "InteractionID (countTagsInRegions) res="+res+" filterThreshold="+threshold+"\tChr\tStart\tEnd\tStrand\tChr\tStart\tEnd";
		for(String dir: dirNames)
			header+="\t"+dir;
		System.out.println(header);
		for(String chrs: map.keySet())
		{
			String[] csplit = chrs.split(",");
			for(String positions : map.get(chrs).keySet())
			{
				double min = 9999999999999.0;
				double max = 0;
				int nonzeros = 0;
				for(i = 0; i < dirs; i++)
				{
					double val = (1000000.0*map.get(chrs).get(positions)[i])/totalCounts[i];
					if(val > 0) nonzeros++;
					if(val > max) max = val;
					if(val < min) min = val;
				}
				boolean badRepValues = false;
				if(expNames.size() > 0)
				{
					//lets go through and check
					for(String exp: expNames.keySet())
					{
						double emax = 0;
						double emin = 99999999;
						for(int index: expNames.get(exp))
						{
							double val = (1000000.0*map.get(chrs).get(positions)[index])/totalCounts[index];
							if(val > emax) emax = val;
							if(val < emin) emin = val;
						}
						if((emax-emin)/(emax+1.0+emin) > maxRepThresh)
						{
							badRepValues=true;
							break;
						}
							
					}
				}
				if (nonzeros > 0 && !badRepValues && (max-min)/(max+1.0+min) > 0.1)
				{
					String[] psplit = positions.split(",");
					int start = Integer.parseInt(psplit[0])*res;
					int end = Integer.parseInt(psplit[1])*res;
					String out = "Interaction"+count+"\t"+csplit[0]+"\t"+start+"\t"+(start+res)+"\t+\t"+csplit[1]+"\t"+end+"\t"+(end+res);
					for(i = 0; i < dirs; i++)
						out+="\t"+map.get(chrs).get(positions)[i];
					System.out.println(out);
					count++;
				}
			}
		}
		System.err.println(count);
	}
	
	private int loadIn(File f, int nameIndex) throws FileNotFoundException
	{
		SuperScanner s = new SuperScanner(f.getAbsolutePath());
		int count = 0;
		while(s.hasMore())
		{
			String line = s.getLine();
			String[] split = line.split("\t");
			String chr1 = split[1];
			String chr2 = split[6];
			int start = Integer.parseInt(split[2])/res;
			int end = Integer.parseInt(split[7])/res;
			if(split[1].compareTo(split[6]) > 0)
			{
				chr1 = new String(split[6]);
				chr2 = new String(split[1]);
				int temp = end;
				end = start;
				start = temp;
			}
			if(chr1.compareTo(chr2) == 0 && start==end)
				continue;  //skip this tag because it is a self-ligation
			//otherwise lets add it
			if(map.get(""+chr1+","+chr2) == null) map.put(""+chr1+","+chr2,new HashMap<String,int[]>());
			if(map.get(""+chr1+","+chr2).get(start+","+end) == null) map.get(""+chr1+","+chr2).put(start+","+end, new int[dirs]);
			map.get(""+chr1+","+chr2).get(start+","+end)[nameIndex]++;
			if(count++ % 1000000 == 0)
				System.err.print(".");
		}
		
		s.close();
		return count;
	}

	public static void main(String[] args) {
		new countTagsInRegions(args);

	}

}
