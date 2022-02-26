import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Scanner;


public class grabFasta {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 2)
		{
			System.err.println("Usage: grabFasta fasta.fa names.txt > filtered.fa  Works with GZIPPED input");
			System.exit(1);
		}
		HashSet<String> names = new HashSet<String>();
		SuperScanner s;
		s = new SuperScanner(args[1]);

		while(s.hasMore())
		{
			String line = s.getLine();
			names.add(line);
		}
		s.close();
		s = new SuperScanner(args[0]);
		String line = "";
		boolean good = false;
		int count = 0; 
		while(s.hasMore())
		{
			line = s.getLine();
			if(line.startsWith(">"))
			{
				String[] ids = line.substring(1).split(" ");
				boolean found = false;
				for(String id: ids)
				{
					if(names.contains(id))
					{
						good = true;
						System.out.println(line);
						found = true;
						count++;
						break;
					}
				}
				if(!found) good = false;
			}
			else if (good)
			{
				System.out.println(line);
			}
		}
		System.err.println("Found "+count+" out of "+names.size());
		
	}

}
