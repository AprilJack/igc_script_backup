import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


public class GrabGTF {

	public GrabGTF(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: GrabGTF fromGTF.gtf list.txt > out.gtf\n"
					+ "Empty lines are disregarded");
			System.exit(1);
		}
		HashSet<String> lines = new HashSet<String>();
		SuperScanner s = new SuperScanner(args[1]);
		while(s.hasMore())
		{
			String line = s.getLine();
			if(line.length() > 0)
				lines.add(line);
		}
		s.close();
		s = new SuperScanner(args[0]);
		int count = 0;
		int ecount = 0;
		while(s.hasMore())
		{
			String line = s.getLine();
			if(line.contains("transcript_id"))
			{
				String id = line.substring(line.indexOf("transcript_id")+15,line.length());
				id = id.substring(0,id.indexOf('"'));
				if(lines.contains(id))
				{
					System.out.println(line);
					ecount++;
				}
			}
			count++;
			if(count%100000==0)
				System.err.print(".");
		}
		System.err.println(" had "+ecount);
		s.close();
		
	}
	
	private LinkedList<String> getSpectra(String txt, int size)
	{
		LinkedList<String> result = new LinkedList<String>();
		for(int i = 0; i < txt.length()-size; i++)
		{
			result.add(txt.substring(i, i+size));
		}
		return result;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new GrabGTF(args);
	}

}
