import java.util.HashMap;
import java.util.LinkedList;


public class GTF {
	HashMap<String,LinkedList<Annotation>> byChr = new HashMap<String,LinkedList<Annotation>>();
	HashMap<String,Annotation> byId = new HashMap<String,Annotation>();
	
	/**
	 * Sorts the byChr Annotations to ensure that they are ordered by their left-most coord
	 */
	public void sort()
	{
		
	}
	
	public GTF(String gtfFile)
	{
		SuperScanner ss = new SuperScanner(gtfFile);
		System.err.print("Loading "+gtfFile);
		int lineCount = 0;

		while(ss.hasMore())
		{
			lineCount++;
			String line = ss.getLine();
			String[] split = line.split("\t");
			if(split.length > 8 && split[2].length() >0)
			{
				String chrStr = split[6]+split[0];
				String ann = split[2];
				if(ann.compareTo("exon")==0)
				{
					int start = Integer.parseInt(split[3].trim());
					int end = Integer.parseInt(split[4].trim());
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					if(byId.get(id) == null)
					{
						Annotation an = new Annotation(chrStr,id);
						an.addRegion(start, end);
						if(byChr.get(chrStr)==null)
							byChr.put(chrStr,new LinkedList<Annotation>());
						byChr.get(chrStr).add(an);
						byId.put(id, an);
						an.addAnnotation(split[8]);
					}
					else
					{
						byId.get(id).addRegion(start,end);
					}
				}
			}
			if(lineCount%10000 == 0) System.err.print(".");
		}
		System.err.println();
	}

}
