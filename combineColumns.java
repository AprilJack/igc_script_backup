import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;


public class combineColumns {

	public static void main(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Concatenates the specified columns from the specified files.");
			System.err.println("Usage: combineColumns file1.csv file2.tab ... range");
			System.err.println("Range can include , or - as well as just a value: 1,3-4,6");
			System.exit(1);
		}
		// TODO Auto-generated method stub
		LinkedList<LinkedList<String>> columns = new LinkedList<LinkedList<String>>();
		String which = args[args.length-1];
		HashSet<Integer> range = new HashSet<Integer>();
		processRanges(which,range);
		for(int i = 0; i < args.length-1; i++)
		{
			LinkedList<LinkedList<String>> fcolumns = new LinkedList<LinkedList<String>>();
			try {
				Scanner s = new Scanner(new File(args[i]));
				for(int c: range)
				{
					 //init
					fcolumns.add(new LinkedList<String>());
				}
				while(s.hasNextLine())
				{
					String line = s.nextLine();
					String[] split = null;
					if(line.indexOf('\t') != -1)
						split= line.split("\t");
					else 
						split= line.split(",");
					int j = 0;
					for(int c: range)
					{
						try{
							fcolumns.get(j).add(split[c-1]);
							j++;
						}
						catch(Exception e)
						{fcolumns.get(j).add("0");}
					}
				}
				columns.addAll(fcolumns);
				s.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		int row = 0;
		while(row < columns.get(0).size())
		{
			for(int i = 0; i < columns.size()-1; i++)
			{
				System.out.print(columns.get(i).get(row)+"\t");
			}
			System.out.println(columns.get(columns.size()-1).get(row));
			row++;
		}
	}
	
	public static void processRanges(String range, HashSet<Integer> list)
	{
		if(range.indexOf('-')!= -1)
		{
			String[] split = range.split("-");
			int min = Integer.parseInt(split[0]);
			int max = Integer.parseInt(split[1]);
			for(int i = min; i <= max; i++)
				list.add(i);
		}
		else if(range.indexOf(',') != -1)
		{
			String[] vals = range.split(",");
			for(String s: vals)
				processRanges(s,list);
		}
		else
			list.add(Integer.parseInt(range));
	}

}
