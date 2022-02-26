import java.io.File;
import java.util.LinkedList;


public class removePeakTags {
	
	public removePeakTags(String[] args)
	{
		if(args.length < 1)
		{
			System.err.println("Usage: removePeakTags tags.tsv [size=1000] [threshold=100] > tags_filtered.tsv \n"
					+ "Filters Homer tag files. If threshold count tags occur within size bases, we assume a peak and stop outputing.");
			System.exit(1);
		}
		int size = 1000;
		int threshold = 100;
		if(args.length > 1)
			size = Integer.parseInt(args[1]);
		if(args.length > 2)
			threshold = Integer.parseInt(args[2]);
		SuperScanner ss = new SuperScanner(args[0]);
		LinkedList<String[]> plus = new LinkedList<String[]>();
		LinkedList<String[]> minus = new LinkedList<String[]>();
		boolean plusPeak = false;
		boolean minusPeak = false;
		while(ss.hasMore())
		{
			String line = ss.getLine();
			String[] split = line.split("\t");
			if(split[3].compareTo("1")==0)
			{
				minus.add(split);
				if(minus.size() > threshold)
				{
					//we must check to see if want to output
					int start = Integer.parseInt(minus.getFirst()[2]);
					int last = Integer.parseInt(minus.getLast()[2]);
					if(last-start < size)
					{
						if(!minusPeak)
							minusPeak = true;
						else
							minus.removeFirst();  //no hope
					}
					else
					{
						if(minusPeak)
						{
							//remove all but the last. 
							minus.clear();
							minus.add(split);
							minusPeak = false;
						}
						String[] first = minus.removeFirst();
						String compiled = first[0];
						for(int i = 1; i < first.length; i++)
							compiled+="\t"+first[i];
						System.out.println(compiled);
					}
				}
			}
			else
			{
				plus.add(split);
				if(plus.size() > threshold)
				{
					//we must check to see if want to output
					int start = Integer.parseInt(plus.getFirst()[2]);
					int last = Integer.parseInt(plus.getLast()[2]);
					if(last-start < size)
					{
						if(!plusPeak)
							plusPeak = true;
						else
							plus.removeFirst();//no chance for that line
					}
					else
					{
						if(plusPeak)
						{
							//remove all but the last. 
							plus.clear();
							plus.add(split);
							plusPeak = false;
						}
						String[] first = plus.removeFirst();
						String compiled = first[0];
						for(int i = 1; i < first.length; i++)
							compiled+="\t"+first[i];
						System.out.println(compiled);
					}
				}
			}

		}
		if(!plusPeak)
		{
			for(String[] first: plus)
			{
				String compiled = first[0];
				for(int i = 1; i < first.length; i++)
					compiled+="\t"+first[i];
				System.out.println(compiled);
			}
		}
		if(!minusPeak)
		{
			for(String[] first: minus)
			{
				String compiled = first[0];
				for(int i = 1; i < first.length; i++)
					compiled+="\t"+first[i];
				System.out.println(compiled);
			}
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new removePeakTags(args);
	}

}
