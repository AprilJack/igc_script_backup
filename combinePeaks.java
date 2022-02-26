import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;


public class combinePeaks {

	HashMap<String,LinkedList<Peak>> peaks = new HashMap<String,LinkedList<Peak>>();

	boolean stranded = false;
	LinkedList<String> theFiles = new LinkedList<String>();
	double[] mins = null;
	
	public combinePeaks(String[] files) {
		if(files.length < 1)
		{
			System.err.println("Usage: combinePeaks [-stranded] [-countCol=N] peaks1.txt peaks2.txt ...");
			System.err.println("Will combine findPeak peaks that overlap in the files and output their respective normalized counts to std out");
			System.err.println("First arg should be -stranded to specify that you want to keep peaks on separate strands separate\n"
					+ "You can also use -countCol=N to specify a different 0-based column for counts (eg. column 8 for raw tag counts)");
			System.exit(1);
		}
		mins = new double[files.length];
		for(int i = 0; i < mins.length; i++)
		{
			mins[i]=Double.POSITIVE_INFINITY;
		}
		int index = 0;
		int countCol = 5;
		for(String file: files)
		{
			if(file.compareTo("-stranded")==0)
			{
				stranded=true;
				System.err.println("Assuming stranded peaks!");
				continue;
			}
			if(file.startsWith("-countCol="))
			{
				countCol = Integer.parseInt(file.substring(10));
				System.err.println("Will use column: "+countCol);
				continue;
			}
			theFiles.add(file);
			SuperScanner fs = new SuperScanner(file);
			int newAdded = 0;
			while(fs.hasMore())
			{
				String line = fs.getLine();
				if(!line.startsWith("#"))
				{
					String[] split = line.split("\t");
					if(split.length > 5)
					{
						int start = Integer.parseInt(split[2]);
						int end = Integer.parseInt(split[3]);
						double count = Double.parseDouble(split[countCol]);
						if(count < mins[index])
							mins[index]=count;
						boolean plus = true;
						if(split[4].compareTo("-")==0)
							plus = false;
						boolean added = false;
						if(peaks.get(split[1])!= null)
						{
							for(Peak p: peaks.get(split[1]))
							{
								if(p.addPeak(split[1], plus|!stranded, start, end, file, count))
								{
									added =true;
									break;
								}
							}
						}
						if(!added)
						{
							newAdded++;
							Peak p = new Peak(split[1],plus|!stranded,start,end,file,count);
							if(peaks.get(split[1])== null) peaks.put(split[1], new LinkedList<Peak>());
							peaks.get(split[1]).add(p);
						}
					}
				}
			}
			System.err.println("Parsed "+file+" containing "+newAdded+" new peaks!");
			index++;
		}
		//now lets output the merged peak counts
		StringBuilder sb = new StringBuilder();
		LinkedList<String> chrs = new LinkedList<String>();
		for(String chr: peaks.keySet())
			chrs.add(chr);
		Collections.sort(chrs);
		int totalPeaks = 0;
		System.out.print("#Peak\tchr\tstart\tend\tstrand");
		for(String f: theFiles)
			System.out.print("\t"+f);
		System.out.println();
		for(String chr: chrs)
		{
			int thusFar = 1;
			Collections.sort(peaks.get(chr));
			for(Peak p: peaks.get(chr))
			{
				totalPeaks++;
				sb.append(p.chr+"_"+thusFar+"\t"+p.chr+"\t"+p.start+"\t"+p.end+"\t"+p.getStrand());
				int i = 0;
				for(String file: theFiles)
				{
					if(p.counts.get(file) != null)
						sb.append("\t"+p.counts.get(file));
					else sb.append("\t"+mins[i]);
					i++;
				}
				sb.append(System.lineSeparator());
				thusFar++;
			}
		}
		System.out.println(sb.toString());
		System.err.println("Finished combining "+totalPeaks+" peaks!");
	}
	
	private class Peak implements Comparable<Peak>
	{
		String chr;
		boolean plus;
		int start;
		int end;
		HashMap<String,Double> counts = new HashMap<String,Double>();
		
		public Peak(String chr, boolean plus,int start2, int end2, String file,
				double count) {
			this.chr = chr; this.start = start2; this.end = end2; counts.put(file,count);
			this.plus = plus;
		}

		public String getStrand() {
			if(plus) return "+";
			else return "-";
		}

		public boolean addPeak(String chr, boolean plus, int start, int end, String file, double count)
		{
			
			if(this.chr.compareTo(chr) == 0 && this.plus==plus)
			{
				int size = (end-start)+(this.end-this.start); //lets extend the peak we are testing this one against
				if((this.start >= start-size && this.start <= end+size)||(this.end >= start-size && this.end <= end+size)||
						(start >= this.start-size && start <= this.end+size)||(end >= this.start-size && end <= this.end+size))
				{
					//peaks overlap
					this.start = Math.min(start,this.start);
					this.end = Math.max(end,this.end);
					counts.put(file,count);
					return true;
				}
			}
			return false;
		}


		@Override
		public int compareTo(Peak p) {
			// TODO Auto-generated method stub
			if(p.start < start)
				return -1;
			else if(p.start > start)
				return 1;
			else return 0;
		}
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new combinePeaks(args);
	}

}
