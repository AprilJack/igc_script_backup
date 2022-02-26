import java.util.HashMap;


public class getUniqueORFs {

	public getUniqueORFs(String[] args) {
		if(args.length == 0)
		{
			System.err.println("Will filter out only the unique ORFs in an annotated list (TranscriptAnnotater). \n"
					+ "ORFs with the same chr and coords will be collapsed into one and output to stdout. Searches for Chr, Start, and End columns");
			System.exit(1);
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
		}
		int count = 0;
		int unique = 0;
		while(s.hasMore())
		{
			String line = s.getLine();
			String[] split = line.split("\t");
			String id = split[chrCol]+","+split[startCol]+","+split[endCol];
			if(map.get(id) == null){
				map.put(id,line);
				System.out.println(line);
				unique++;
			}
			count++;
			if(count%100000 ==0) System.err.print(".");
		}
		System.err.println(" Total: "+count+" Unique: "+unique);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new getUniqueORFs(args);
	}

}
