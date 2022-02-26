import java.awt.Toolkit;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


public class sortPairedFastqs {

	public sortPairedFastqs(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: sortPairedFastqs fastq_R1.fastq.gz fastq_R2.fastq.gz");
			System.err.println("Outputs fastq_sorted_R1.fastq.gz fastq_sorted_R2.fastq.gz");
			System.exit(1);
		}
		try {
			Scanner ss = new Scanner(new GZIPInputStream(new FileInputStream(args[0])));
			Scanner ss2 = new Scanner(new GZIPInputStream(new FileInputStream(args[1])));
			PrintWriter pw1 = new PrintWriter(new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(args[0].substring(0,args[0].length()-9)+"_fixed.fastq.gz"))));
			PrintWriter pw2 = new PrintWriter(new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(args[1].substring(0,args[1].length()-9)+"_fixed.fastq.gz"))));
			System.err.println("Outputing matching reads...");
			HashMap<String,String[]> map = new HashMap<String,String[]>();
			int count = 0;
			while(ss.hasNextLine() && ss2.hasNextLine())
			{
				String[] result = new String[4];
				for(int i = 0; i < result.length; i++)
				{
					String line = ss.nextLine();
					result[i]=line;
					//System.err.println(line);
				}
				String[] result2 = new String[4];
				for(int i = 0; i < result2.length; i++)
					result2[i]=ss2.nextLine();
				if(map.get(result2[0].split(" ")[0])!= null)
				{
					String[] stored = map.get(result2[0].split(" ")[0]);
					for(int i = 0; i < stored.length; i++)
						pw1.println(stored[i]);
					for(int i = 0; i < result2.length; i++)
						pw2.println(result2[i]);
					map.remove(result2[0].split(" ")[0]);
				}
				else 
					map.put(result2[0].split(" ")[0],result2);
				if(map.get(result[0].split(" ")[0])!= null)
				{
					String[] stored = map.get(result[0].split(" ")[0]);
					for(int i = 0; i < stored.length; i++)
						pw2.println(stored[i]);
					for(int i = 0; i < result.length; i++)
						pw1.println(result[i]);
					map.remove(result[0].split(" ")[0]);
				}
				else 
					map.put(result[0].split(" ")[0],result);
				count++;
				if(count %1000000 == 0)
					System.err.println((count)/1000000+" M reads. Used:  "+(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024)+" Mb");
				
			}
			System.err.println();
			ss.close();
			pw1.close();
			pw2.close();
			ss2.close();		
			System.err.println();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

	public static void main(String[] args) {
		new sortPairedFastqs(args);
	}

}
