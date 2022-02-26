import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;


public class processMergedInteractions {

	public processMergedInteractions(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Takes a basic interactions file and the original mergedInteractions file. \n"
					+ "Outputs a peak file containing end points for the interactions.\n"
					+ "\nUsage: processMergedInteractions interactions.txt mergedInteractions.txt > endpoint_peaks.txt");
			System.exit(1);
		}
		try {
			Scanner s = new Scanner(new File(args[1]));
			HashMap<String,String[]> map = new HashMap<String,String[]>();
			s.nextLine(); //get rid of header for the merged file 
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\t");
				map.put(split[0],split);
			}
			s.close();
			s = new Scanner(new File(args[0]));
			String[] hsplit = s.nextLine().split("\t"); //get rid of header for the merged file 
			String header ="PeakID (processMergedInteractions "+args[0]+" "+args[1]+")\tChr\tStart\tEnd\tStrand";
			for(int i = 1; i < hsplit.length; i++)
				header+="\t"+hsplit[i];
			System.out.println(header);
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\t");
				if(map.get(split[0])!= null)
				{
					//lets append
					String out = split[0]+"_1\t"+map.get(split[0])[1]+"\t"+map.get(split[0])[2]+"\t"+map.get(split[0])[3]+"\t"+map.get(split[0])[4];
					for(int i = 1; i < split.length; i++)
						out+="\t"+split[i]; //add the normalized counts and p-values back at the end
					System.out.println(out);
					out = split[0]+"_2\t"+map.get(split[0])[5]+"\t"+map.get(split[0])[6]+"\t"+map.get(split[0])[7]+"\t"+map.get(split[0])[4];
					for(int i = 1; i < split.length; i++)
						out+="\t"+split[i]; //add the normalized counts and p-values back at the end
					System.out.println(out);
				}
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new processMergedInteractions(args);
	}

}
