import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;


public class fastqDictionary {

	HashMap<String,HashMap<String,Integer>> map = new HashMap<String,HashMap<String,Integer>>();
	
	public fastqDictionary(String[] args) throws IOException
	{
		if(args.length < 1)
		{
			System.err.println("Usage: fastqDictionary read.fastq[.gz] > dictionary.fa"
					+ "\nReads in the fastq file and extracts the barcode (seq in name after last :). Then counts up the occurances of the barcode with a specific read prefix."
					+ "\nAssumes you already ran extractUMI and corrected for errors in the read (e.g. with karect)"
					+ "\nOutputs a table in tsv format to std out");
			System.exit(1);
		}
		SuperScanner ss = new SuperScanner(args[0]);
		int count = 0;
		while(ss.hasMore())
		{
			String name = ss.getLine();
			String barcode = name.split(":")[name.split(":").length-1];
			String seq = ss.getLine();
			String prefix = seq.substring(0,Math.min(seq.length(), 30));
			if(map.get(barcode) == null)
				map.put(barcode, new HashMap<String,Integer>());
			if(map.get(barcode).get(prefix)== null)
				map.get(barcode).put(prefix, 1);
			else
				map.get(barcode).put(prefix, map.get(barcode).get(prefix));
			ss.getLine(); //note
			ss.getLine(); //quality (for now forget about it)
			count++;
			if(count % 400000 == 0)
				System.err.print(".");
		}
		System.err.println("Done");
		System.out.println("Barcode\tTop Seq\tTop Seq #\tSecond Seq\tSecond Seq #");
		for(String barcode: map.keySet())
		{
			HashMap<String,Integer> hits = map.get(barcode);
			if(hits != null)
			{
				int best=0;
				int second=0;
				String bestHit="";
				String secondHit="";
				for(String hit: hits.keySet())
				{
					if(hits.get(hit)> best)
					{
						second=best;
						secondHit=bestHit;
						best=hits.get(hit);
						bestHit=hit;
					}
				}
				System.out.println(barcode+"\t"+bestHit+"\t"+best+"\t"+secondHit+"\t"+second);
			}
			else System.err.println("Should not happen!");
		}
		ss.close();
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new fastqDictionary(args);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
