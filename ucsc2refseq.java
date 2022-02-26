

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

public class ucsc2refseq {

	public static void main(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage:ucsc2refseq fileWithUCSCids ucsc2refseq.conversion");
			System.err.println("Only keeps the lines that have valid mappings. I.e. if the second column is missing in the coversion file");
			System.exit(1);
		}
		// TODO Auto-generated method stub
		//load in a ucsc2refseq file and a file containing ucsc gene ids 
		//keep the top expressing one
		HashMap<String,String> map = new HashMap<String,String>(); //ucsc->refseq
		
		try {
			Scanner s = new Scanner(new File(args[1]));
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\t"); //split on tab
				if(split.length > 1)
					map.put(split[0], split[1]);
				else
					map.put(split[0],"!");
			}
			s.close();
			s = new Scanner(new File(args[0]));
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("[\t\";]+");
				boolean print = true;
				for(String str: split)
				{
					if(map.get(str) != null)
					{
						if(map.get(str).compareTo("!")==0)
							print = false;
						line = line.replaceFirst(str, map.get(str));
					}
				}
				if(print)
					System.out.println(line);
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}

}
