import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.zip.GZIPOutputStream;


public class FastqDownsampler {

	public FastqDownsampler(String[] args) throws IOException
	{
		if(args.length < 2)
		{
			System.err.println("Usage: FastaDownsampler read1.fastq[.gz] [read2.fastq[.gz]] prob_to_keep"
					+ "\nGenerates read1_prob.fastq.gz and read2_prob.fastq.gz in same location");
			System.exit(1);
		}
		//lets figure out the prob first
		double prob = Double.parseDouble(args[args.length-1]);
		SuperScanner ss = new SuperScanner(args[0]);
		SuperScanner ss2 = null;
		BufferedWriter pw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args[0].replace(".fastq", "_"+prob+".fastq")))));
		BufferedWriter pw2= null;
		if(args.length > 2)
		{
			ss2 = new SuperScanner(args[1]);
			pw2 =new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args[1].replace(".fastq", "_"+prob+".fastq")))));
		}
		//we assume that if there are paired reads then they are in the same order and that there is an equal number
		int count = 0;
		while(ss.hasMore())
		{
			//lets figure out if we want to keep this read
			if(Math.random() < prob)
			{
				pw.append(ss.getLine()); pw.append(System.lineSeparator());
				pw.append(ss.getLine()); pw.append(System.lineSeparator());
				pw.append(ss.getLine()); pw.append(System.lineSeparator());
				pw.append(ss.getLine()); pw.append(System.lineSeparator());
				if(pw2 != null)
				{
					pw2.append(ss2.getLine()); pw2.append(System.lineSeparator());
					pw2.append(ss2.getLine()); pw2.append(System.lineSeparator());
					pw2.append(ss2.getLine()); pw2.append(System.lineSeparator());
					pw2.append(ss2.getLine()); pw2.append(System.lineSeparator());
				}
			}
			else
			{
				ss.getLine();ss.getLine();ss.getLine();ss.getLine();
				if(ss2 != null)ss2.getLine();ss2.getLine();ss2.getLine();ss2.getLine();
			}
			count++;
			if(count % 400000 == 0)
				System.err.print(".");
		}
		System.err.println("Done");
		ss.close();
		pw.close();
		if(ss2 != null)
			ss2.close();
		if(pw2!= null)
			pw2.close();
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new FastqDownsampler(args);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
