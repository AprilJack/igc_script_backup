import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;


public class fastqSplitter {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length == 0)
		{
			System.err.println("Usage: fastqSplitter paired_sequential_reads.fastq Prefix");
			System.err.println("Outputs in gzipped format but allows for both gzipped and unzipped input");
			System.exit(1);
		}
		SuperScanner ss = new SuperScanner(args[0]);
		try {
			BufferedOutputStream out = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(args[1]+"_R1.fastq.gz")),8096);
			BufferedOutputStream err = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(args[1]+"_R2.fastq.gz")),8096);

			while(ss.hasMore())
			{
				out.write(ss.getLine().getBytes());out.write(System.lineSeparator().getBytes());
				out.write(ss.getLine().getBytes());out.write(System.lineSeparator().getBytes());
				out.write(ss.getLine().getBytes());out.write(System.lineSeparator().getBytes());
				out.write(ss.getLine().getBytes());out.write(System.lineSeparator().getBytes());
				err.write(ss.getLine().getBytes());err.write(System.lineSeparator().getBytes());
				err.write(ss.getLine().getBytes());err.write(System.lineSeparator().getBytes());
				err.write(ss.getLine().getBytes());err.write(System.lineSeparator().getBytes());
				err.write(ss.getLine().getBytes());err.write(System.lineSeparator().getBytes());
			}
			ss.close();
			out.close();
			err.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
