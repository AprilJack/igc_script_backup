import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


public class Grab {

	public Grab(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: Grab list.txt fromThisFile.file > out.txt");
			System.exit(1);
		}
		HashSet<String> lines = new HashSet<String>();
		SuperScanner s = new SuperScanner(args[0]);
		while(s.hasMore())
		{
			String line = s.getLine();
			lines.add(line);
		}
		s.close();
		s = new SuperScanner(args[1]);
		int count = 0;
		while(s.hasMore())
		{
			String line = s.getLine();
			if(line.contains("transcript_id"))
			{
				String id = line.substring(line.indexOf("transcript_id")+15,line.length());
				id = id.substring(0,id.indexOf('"'));
				if(lines.contains(id))
					System.out.println(line);
			}
			count++;
			if(count%100000==0)
				System.err.print(".");
		}
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
		new Grab(args);
	}

}
