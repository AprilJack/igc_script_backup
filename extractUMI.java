import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.zip.GZIPOutputStream;


public class extractUMI {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//we are reading in a file one line at a time
		if(args.length < 2)
		{
			System.err.println("Usage: saveUMI N fastq.gz > fastq.gz");
			System.err.println("Cuts and inserts the first N bases of the read into its name.\n"
					+ "Output to standard out as gzipped");
			System.exit(1);
		}
		int N = Integer.parseInt(args[0]);
		try {
			SuperScanner s = new SuperScanner(args[1]);
			PrintWriter pw = new PrintWriter(new GZIPOutputStream(System.out));
			int lines = 0;
			StringBuilder sb = new StringBuilder(1024*1024*100);
			int totalReads = 0;
			while(s.hasMore())
			{
				totalReads++;
				if(lines++ % 100000 == 0)
				{
					pw.print(sb.toString());
					sb = new StringBuilder(1024*1024*100);
					System.err.println("Processed "+(lines)+" reads in file "+args[1]+"...");
				}

				String name = s.getLine();
				String seq = s.getLine();
				String plus = s.getLine(); // + string
				String qual = s.getLine();
				//extract UMI
				String UMI = seq.substring(0,N);
				String UMIqual = qual.substring(0,N);
				seq = seq.substring(N);
				qual = qual.substring(N);
				sb.append(name); sb.append(UMI);
				//sb.append(line.substring(1));
				sb.append(System.lineSeparator());
				sb.append(seq);sb.append(System.lineSeparator());
				sb.append(plus);sb.append(System.lineSeparator());
				sb.append(qual);sb.append(System.lineSeparator());
			}
			pw.print(sb.toString());
			s.close();
			pw.close();System.err.println("DONE! Processed "+totalReads);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
