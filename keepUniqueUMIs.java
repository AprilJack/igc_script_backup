import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;


public class keepUniqueUMIs {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 1)
		{
			System.err.println("Keeps only the unique UMIs for each position (Assume name of read is just the 5letter UMI). Outputs sam file");
			System.err.println("Usage: keepUniqueUMIs samFile1 [samFile2] [samFile3] ...");
			System.exit(1);
		}
		for(String s: args)
		{
			if(s.endsWith(".sam"))
			{
				processFile(s);
			}
		}
	}

	static final int capacity = 100*1024*1024;
	private static void processFile(String str) {
		// TODO Auto-generated method stub
		Scanner s;
		HashMap<String,HashMap<String,Short>> hs =new HashMap<String,HashMap<String,Short>>();
		try {
			PrintWriter pw = new PrintWriter(str.substring(0,str.indexOf(".sam"))+"_uniqueUMI.sam");
			PrintWriter pw2 = new PrintWriter(str.substring(0,str.indexOf(".sam"))+"_summary.csv");
			s = new Scanner(new File(str));
			StringBuilder sb = new StringBuilder(capacity);
			int count =0;
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				if(line.startsWith("@"))
				{
					sb.append(line);sb.append(System.lineSeparator());
					continue;
				}
				String[] split = line.split("\t");
				int pos = Integer.parseInt(split[3]);
				String flag = Integer.toBinaryString(Integer.parseInt(split[1]));
				boolean reverse = false;
				if(flag.length() > 4 && flag.charAt(flag.length()-5)=='1')
					reverse = true;
				if(reverse)
					pos+=split[9].length();
				String position = split[2]+","+pos;	
				String umi = split[0];
				if(hs.get(position) == null)
				{
					hs.put(position, new HashMap<String,Short>());
				}
				if(hs.get(position).get(umi) == null)
				{
					hs.get(position).put(umi, (short)1);	
					if(sb.capacity() < line.length()+2)
					{
						System.err.println(str+": processed "+(count++));
						pw.print(sb.toString());
						sb = new StringBuilder(capacity);
					}
					sb.append(line);sb.append(System.lineSeparator());
				}
				else
				{
					hs.get(position).put(umi,(short)( hs.get(position).get(umi)+1));
				}

			}
			pw.print(sb.toString());
			s.close();
			pw.close();
			int[] UMIhisto = new int[1024];
			for(String key: hs.keySet())
			{
				int cc = hs.get(key).size();
				if(cc < 1025)
					UMIhisto[cc-1]++;
				else
					System.err.println("Encountered more than 1024 unique UMIs for a position!");
			}
			pw2.println("#UniqueUMIs,Positions");
			for(int i = 0; i < UMIhisto.length; i++)
			{
				pw2.println((i+1)+"\t"+UMIhisto[i]);
			}
			pw2.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println(str+"\t"+hs.size());
		
	}

}
