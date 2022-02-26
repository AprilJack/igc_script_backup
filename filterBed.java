
import java.io.FileNotFoundException;
import java.io.PrintWriter;


public class filterBed {

	public static void main(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: filterBed length bedfile1 [bedfile2]...");
			System.err.println("Outputs to new file with \"_filtered.bed\" suffix");
			System.exit(1);
		}
		int threshold = Integer.parseInt(args[0]);
		for(int i = 1; i < args.length; i++)
		{
			FastScanner fs = new FastScanner(args[i]);
			try {
				PrintWriter pw = new PrintWriter(args[i].substring(0, args[i].indexOf("."))+"_filtered.bed");
				while(fs.hasMore())
				{
					String line = fs.getLine();
					if(!line.startsWith("#"))
					{
						String[] split = line.split("\t");
						if(split.length > 2)
						{
							long start = Long.parseLong(split[1]);
							long stop = Long.parseLong(split[2]);
							int diff = (int) Math.abs(start-stop);
							if(diff > threshold)
								pw.println(line);
						}
					}
					else pw.println(line);
				}
				pw.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	}
}
