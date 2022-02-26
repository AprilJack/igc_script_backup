import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;


public class orfCombiner {

	HashMap<String,boolean[]> map = new HashMap<String,boolean[]>();
	HashMap<String,String> unique = new HashMap<String,String>();
	public orfCombiner(String[] args) throws FileNotFoundException {
		//assume each arg is a filename that has ORFs in the format outputed by the RunRiboSeq pipeline
		for(int i =0; i < args.length; i++)
		{
			Scanner s = new Scanner(new File(args[i]));
			int count = 0;
			int ucount = 0;
			while(s.hasNextLine())
			{
				count++;
				String line = s.nextLine();
				int tokens = line.split("[+-]").length;
				String loc = line.split("[+-]")[tokens-2]+"|"+line.split("[+-]")[tokens-1].split("_")[0];
				if(line.indexOf('+')==0)
					loc="+"+loc;
				else
					loc="-"+loc;
				if(map.get(loc) == null)
				{
					map.put(loc, new boolean[args.length]);
				}
				if(unique.get(loc)==null)
				{
					unique.put(loc,line);
					ucount++;
				}
				map.get(loc)[i]=true; 
			}
			s.close();
			System.err.println(String.format("Finished reading in %s. Had %d ORFs of which %d were unique thus far.",args[i],count,ucount));
		}
		System.err.println(String.format("Outputing %d unique ORFs",map.keySet().size()));
		//header
		String header = "smORF";
		for(int i = 0; i < args.length; i++)
			header+="\t"+args[i];
		System.out.println(header);
		for(String key: map.keySet())
		{
			String line = unique.get(key);
			for(int i = 0; i < args.length; i++)
			{
				line+="\t"+map.get(key)[i];
			}
			System.out.println(line);
		}
		
	}

	public static void main(String[] args) {
		try {
			if(args.length > 0) {
				new orfCombiner(args);
			} else {
				System.err.println("Usage: orfCombiner fileWithDetectedORFs [fileWithDetectedORFs] ...");
				System.err.println("Will find unique ORFs and output an occupancy matrix.");
			}
				
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
