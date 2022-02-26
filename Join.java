import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;


public class Join {

	public static void main(String[] args) {
		if(args.length == 0)
		{
			System.err.println("Naively joins one or more columns from files. Assumes files have identical number of rows. \n"
					+ "\n\tUsage: combineColumns file1 [file2] ... [fileN] [-sep=,]\n"
					+ "\tuse -sep to replace all instances of the specified sep with a tab.");
		}
		ArrayList<String> lists = new ArrayList<String>();
		String sep = null;
		for(String arg: args)
		{
			if(arg.contains("-sep"))
			{
				sep = arg.substring(5);
				break;
			}
			try {
				Scanner s =new Scanner(new File(arg));
				int i = 0;
				while(s.hasNextLine())
				{
					String line = s.nextLine();
					if(lists.size()== i) lists.add(line);
					else lists.set(i,lists.get(i)+"\t"+line);
					i++;
				}
				s.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		for(String line: lists)
		{
			if(sep == null)
				System.out.println(line);
			else
				System.out.println(line.replaceAll(sep, "\t"));
		}

	}

}
