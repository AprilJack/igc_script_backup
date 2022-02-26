import java.util.HashMap;


public class filterORFs {

	public filterORFs(String[] args) {
		if(args.length == 0)
		{
			System.err.println("Will filter out only the unique ORFs in an annotated list (TranscriptAnnotater). Will also only keep expression values higher than a set threshold.\n"
		
					+ "ORFs with the same chr and coords will be collapsed into one and output to stdout. Searches for Chr, Start, and End columns\n"
					+ "Usage: filterORFs ORFs_ann.txt [col#,Threshold] [col#,Threshold] ...\n"
					+ "Example: filterORFs ORFs.txt 12,0.1 13,1.5 14,0.1 > filtered.txt");
			System.exit(1);
		}
		HashMap<Integer,Double> thresholds = new HashMap<Integer,Double>();
		for(int i = 1; i < args.length; i++)
		{
			String arg = args[i];
			thresholds.put(Integer.parseInt(arg.split(",")[0]),Double.parseDouble(arg.split(",")[1]));
		}
		HashMap<String,String> map = new HashMap<String,String>();
		SuperScanner s = new SuperScanner(args[0]);
		int chrCol,startCol,endCol;
		chrCol=startCol=endCol=0;
		if(s.hasMore())
		{
			String line = s.getLine();
			String[] split = line.split("\t");
			for(int i = 0; i < split.length; i++)
			{
				if(split[i].compareTo("Chr")==0) chrCol =i;
				else if(split[i].compareTo("Start")==0) startCol = i;
				else if(split[i].compareTo("End")==0) endCol = i;
			}
			System.out.println(line);
		}
		int count = 0;
		int unique = 0;
		while(s.hasMore())
		{
			String line = s.getLine();
			String[] split = line.split("\t");
			String id = split[chrCol]+","+split[startCol]+","+split[endCol];
			if(map.get(id) == null){
				boolean passes = true;
				for(int col: thresholds.keySet())
				{
					try{
						double val = Double.parseDouble(split[col]);
						if(val < thresholds.get(col))
						{
							passes = false;
							break;
						}
					}catch(Exception e)
					{
						passes = false;
						break;
					}
				}
				if(passes)
				{
					map.put(id,line);
					System.out.println(line);
					unique++;
				}
			}
			count++;
			if(count%100000 ==0) System.err.print(".");
		}
		System.err.println(" Total: "+count+" Remaining: "+unique);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new filterORFs(args);
	}

}
