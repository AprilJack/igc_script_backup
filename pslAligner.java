import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;


public class pslAligner {

	
	/**
	 * Output pslBig filtered to only contain splLittle entries with the positions added
	 * @param args
	 */
	public pslAligner(String[] args) {
		String pslLittle = args[0];
		String pslBig = args[1];  
		HashMap<String,LinkedList<String>> map = new HashMap<String,LinkedList<String>>();
		try {
			Scanner s = new Scanner(new File(pslLittle));
			if(s.hasNextLine())System.out.println(s.nextLine()); //header
			if(s.hasNextLine())System.out.println(s.nextLine()); //header
			if(s.hasNextLine())System.out.println(s.nextLine()); //header
			if(s.hasNextLine())System.out.println(s.nextLine()); //header
			if(s.hasNextLine())System.out.println(s.nextLine()); //header
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\\s+");
				String target = split[13];
				if(map.get(target) == null) map.put(target, new LinkedList<String>());
				map.get(target).add(line);
			}
			s.close();
			s= new Scanner(new File(pslBig));
			if(s.hasNextLine())s.nextLine(); //header
			if(s.hasNextLine())s.nextLine(); //header
			if(s.hasNextLine())s.nextLine(); //header
			if(s.hasNextLine())s.nextLine(); //header
			if(s.hasNextLine())s.nextLine(); //header
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\\s+");  //split is in pslBIG
				if(map.get(split[9]) != null)
				{
					LinkedList<String> list = map.remove(split[9]);
					for(String a: list)
					{
						String[] split2 = a.split("\\s+");
						
						for(int i = 0; i < split2.length; i++)
						{
							if(i == 0) System.out.print(split2[i]);
							else if(i==8)  //alignment direction to the genome (either + or -)
							{
								if(split2[8].compareTo("++")==0 && split[8].compareTo("+")==0)System.out.print("\t+");
								else if(split2[8].compareTo("++")==0 && split[8].compareTo("-")==0)System.out.print("\t-");
								else if(split2[8].compareTo("-+")==0 && split[8].compareTo("+")==0)System.out.print("\t-");
								else if(split2[8].compareTo("+-")==0 && split[8].compareTo("+")==0)System.out.print("\t-");
								else if(split2[8].compareTo("+-")==0 && split[8].compareTo("-")==0)System.out.print("\t+");
								else if(split2[8].compareTo("-+")==0 && split[8].compareTo("-")==0)System.out.print("\t+");
								else if(split2[8].compareTo("--")==0 && split[8].compareTo("-")==0)System.out.print("\t-");
								else if(split2[8].compareTo("--")==0 && split[8].compareTo("+")==0)System.out.print("\t+");
							}
							else if(i == 13 || i==14)
								System.out.print("\t"+split[i]);
							else if(i == 15)
								System.out.print("\t"+(Long.parseLong(split[i])+Long.parseLong(split2[i])-Long.parseLong(split2[11])*3));
							else if(i == 16)
								System.out.print("\t"+(Long.parseLong(split[15])+Long.parseLong(split2[15])+Long.parseLong(split2[10])*3-Long.parseLong(split2[11])*3));
							else if(i == 17) //qStartsd
								System.out.print("\t1");
							else if(i == 18)
								System.out.print("\t"+(Long.parseLong(split2[10])*3));
							else if(i == 19) //qStarts
								System.out.print("\t0");
							else if(i == 20)
							{
								System.out.print("\t"+(Long.parseLong(split[15])+Long.parseLong(split2[15])-Long.parseLong(split2[11])*3));
							}
							else
								System.out.print("\t"+split2[i]);
						}
						System.out.println();
					}
				}
			}
			for(String str: map.keySet())
			{
				System.err.println(str+" not in the second psl file!");
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args)
	{
		new pslAligner(args);
	}
}
