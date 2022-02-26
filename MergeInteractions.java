import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;


public class MergeInteractions {

	HashMap<String,LinkedList<Interaction>> map = new HashMap<String,LinkedList<Interaction>>();
	HashSet<String> sampleNames = new HashSet<String>();
	int resolution = 0;
	
	public MergeInteractions(String[] args) {
		if(args.length < 2)
		{
			System.err.println("MergeInteractions will merge Homer interaction files such that only the common subset remains. \n"
					+ "Assumes fixed bin boundaries (eg. all bin boundaries are the same between experiments so don't use -center). \n"
					+ "Every other line is ignored since the data is just repeated and flipped. Use interactionsToPeaks to get the endpoints back");
			System.exit(1);
		}
		for(String f: args)
		{
			try {
				Scanner s = new Scanner(new File(f));
				sampleNames.add(f);
				s.nextLine();
				System.err.println("Reading interactions for "+f);
				while(s.hasNextLine())
				{
					String line = s.nextLine();
					String[] split = line.split("\t");
					//lets go through and combine the interactions. 
					String chr1 = split[2];
					String chr2 = split[8];
					int start1 = Integer.parseInt(split[3]);
					int end1 = Integer.parseInt(split[4]);
					resolution = end1-start1;
					int start2 = Integer.parseInt(split[9]);
					int count = (int)Double.parseDouble(split[14]);
					Interaction newi = new Interaction(f,count,chr1,chr2,start1,start2);
					if(chr1.compareTo(chr2) > 0)
					{
						String temp = ""+chr2;
						chr2 = ""+chr1;
						chr1 = temp;
					}
					if(map.get(chr1+chr2) == null) map.put(chr1+chr2, new LinkedList<Interaction>());
					Interaction found = null;
					for(Interaction i: map.get(chr1+chr2))
					{
						if(i.overlapsWith(newi))
						{
							found = i;
							break;
						}
					}
					if(found != null)
					{
						found.integrateInteraction(newi);
					}
					else
						map.get(chr1+chr2).add(newi);
				}
				s.close();

				//now lets output the interactions
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		System.err.println();
		String header ="InteractionID (MergeInteractions)\tChr\tStart\tEnd\tStrand\tChr\tStart\tEnd";
		for(String sample: sampleNames)
			header+="\t"+sample+" tags";
		System.out.println(header);
		int interactionCount = 1;
		for(String key: map.keySet())
		{
			for(Interaction i: map.get(key))
			{
				String out = "Interaction"+interactionCount+"\t"+i.chr1+"\t"+i.start+"\t"+(i.start+resolution)+"\t+\t"+i.chr2+"\t"+i.end+"\t"+(i.end+resolution);
				for(String sample: sampleNames)
				{
					if(i.counts.get(sample) != null)
						out+="\t"+i.counts.get(sample);
					else out+="\t0";
				}
				System.out.println(out);
				interactionCount++;
			}
		}
	}
	
	private class Interaction
	{
		String chr1,chr2;
		int start,end;
		HashMap<String,Integer> counts = new HashMap<String,Integer>();
		
		public Interaction(String id, int count, String chr1, String chr2, int start, int end)
		{
			this.chr1=new String(chr1); this.chr2=new String(chr2); this.start = start; this.end = end;
			counts.put(id, count);
		}
		
		public boolean overlapsWith(Interaction i)
		{
			if(chr1.compareTo(i.chr1)==0 && chr2.compareTo(i.chr2)==0 && start == i.start && end == i.end)
				return true;
			else if(chr1.compareTo(i.chr2)==0 && chr2.compareTo(i.chr1)==0 && start == i.end && end == i.start)
				return true;
			return false;
		}
		

		public void integrateInteraction(Interaction i)
		{
			if(overlapsWith(i))
			{
				for(String key: i.counts.keySet())
				{
					counts.put(key,i.counts.get(key));  //same sample file can only have one count per interaction region
				}
			}
		}
	}
	
	private static final boolean overlaps(int left1,  int right1, int left2,int right2)
	{
		int max = Math.min(right1,right2);
		int min = Math.max(left1,left2);
		return (max-min) > 0; 
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new MergeInteractions(args);
	}

}
