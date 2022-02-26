import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;



public class keepBestpsl {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 2)
		{
			System.err.println("Usage: keepBestPsl fasta_file psl_file\nOutputs two files: fasta_file_processed containing only those seqs which are in the psl, and a psl_processed containing just the best-matched entries");
			System.exit(1);
		}
		FastScanner s = new FastScanner(args[1]);
		try {
			PrintWriter psl = new PrintWriter(args[1]+"_processed.psl");
			PrintWriter fasta = new PrintWriter(args[0]+"_processed.fa");
			HashSet<String> transcripts = new HashSet<String>();
			for(int i = 0; i < 5; i++)
				if(s.max > 0)psl.println(s.getLine()); 
			String current = null;
			int matches = 0;
			String bestMatch = null;
			String line = s.getLine();
			while(s.max > 0)
			{

				String[] split = line.split("\\s+");
				if(current == null || current.compareTo(split[9]) != 0)
				{
					if(current != null)
					{
						psl.println(bestMatch);
						transcripts.add(split[9]);
					}
					current = split[9];
					matches = Integer.parseInt(split[0]);
					bestMatch = line;
				}
				else  //current and split[9] are the same so this could be a better match potentially
				{
					int match_length = Integer.parseInt(split[0]);
					if(match_length > matches)
					{
						bestMatch = line;
						matches = match_length;
					}
				}
				line = s.getLine();
			}
			s.close();
			psl.close();
			s = new FastScanner(args[0]);
			boolean writing = false;
			line = s.getLine();
			while(s.max > 0)
			{
				
				if(line.startsWith(">"))
				{
					writing = false;
					String name = line.substring(1);
					if(!transcripts.contains(name))
					{
						writing = true;
						fasta.println(line);
					}
				}
				else
				{
					if(writing)
						fasta.println(line);
				}
				line = s.getLine();
			}
			s.close();
			fasta.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
