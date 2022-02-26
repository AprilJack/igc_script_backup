

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

public class enrichmentAnalysis {

	public static void main(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: countUpAnnotation <fullHomerAnnotation.txt> <counts_per_element>");
			System.exit(1);
		}
		HashMap<String,Long> tallies = new HashMap<String,Long>();
		HashMap<String,Double> measured = new HashMap<String,Double>();
		HashMap<String,String> replace = new HashMap<String,String>();
		replace.put("N","Intergenic");
		replace.put("E","Exon");
		replace.put("I","Intron");
		replace.put("P","Promoter");
		try {
			Scanner s = new Scanner(new File(args[0]));
			while (s.hasNextLine())
			{
				String[] line = s.nextLine().split("\t");
				if(tallies.get(line[5]) == null) tallies.put(line[5],0l);
				long current = tallies.get(line[5]);
				tallies.put(line[5],current+(Long.parseLong(line[3])-Long.parseLong(line[2])));
			}
			s.close();
			s = new Scanner(new File(args[1]));
			long total = 0;
			for(long l: tallies.values())
				total+=l;
			double measuredTotal = 0;
			while(s.hasNextLine())
			{
				String[] line = s.nextLine().split("\t");
				//lets figure out the expected number of reads for this type of element
				measured.put(line[0],Double.parseDouble(line[1]));
				measuredTotal+=measured.get(line[0]);
			}
			for(String element: tallies.keySet())
			{
				double expected = ((double)tallies.get(element))/total;
				double actual = 0;
				if(measured.get(replace.get(element))!= null)
					element = replace.get(element);
				if(measured.get(element) != null)
					actual = ((double)measured.get(element))/measuredTotal;
				double enrichment = actual/expected;
				if(total == 0 || measuredTotal == 0)
					System.out.println(element+"\t0");
				else
					System.out.println(String.format("%s\t%5.4f",element,enrichment));
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}
