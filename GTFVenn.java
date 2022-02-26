import java.util.HashMap;
import java.util.LinkedList;


public class GTFVenn {

	HashMap<String,HashMap<Integer, LinkedList<Annotation>>> annotations = new HashMap<String,HashMap<Integer, LinkedList<Annotation>>>();
	HashMap<String, Annotation> byId = new HashMap<String,Annotation>();  //assume the ids are all unique across all GTFs
	HashMap<String, int[]> scores = new HashMap<String,int[]>();  //assume the ids are all unique across all GTFs
	HashMap<String, Integer> fileMap = new HashMap<String,Integer>();  //assume the ids are all unique across all GTFs
	public GTFVenn(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: GTFVenn gtf1.txt gtf2.txt [gtf3.txt] ...  > venn.txt\n"
					+ "For each transcript across all gtfs, the score across all other GTFs is recorded.  \n"
					+ "Then transcripts with a score of at leaast 1 (some overlap) are considered a hit \n"
					+ "and a venn diagram of common hits is constructed for visualization.");
			System.exit(1);
		}
		
		for(int gtf = 0; gtf < args.length; gtf++)
		{
			SuperScanner ss = new SuperScanner(args[gtf]);
			System.err.print("Loading "+args[gtf]);
			int lineCount = 0;
			int count = 0;
			while(ss.hasMore())
			{
				lineCount++;
				String line = ss.getLine();
				String[] split = line.split("\t");
				if(split.length > 8 && split[2].length() >0)
				{
					String chrStr = split[0];
					String ann = split[2];
					if(ann.compareTo("exon")==0)
					{
						int start = Integer.parseInt(split[3].trim());
						int end = Integer.parseInt(split[4].trim());
						String id = split[8].replaceAll("\"", "");
						id=id.replaceAll("\"", "");
						id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
						String key = chrStr;
						if(byId.get(id) == null)
						{
							byId.put(id,new Annotation(chrStr,id));
							scores.put(id, new int[args.length]);
							fileMap.put(id,gtf);
							if(annotations.get(chrStr)==null) annotations.put(chrStr,new HashMap<Integer,LinkedList<Annotation>>());
							if(annotations.get(chrStr).get(gtf) == null) annotations.get(chrStr).put(gtf, new LinkedList<Annotation>());
							annotations.get(chrStr).get(gtf).add(byId.get(id));
							count++;
							
						}
						byId.get(id).addRegion(start, end);
						byId.get(id).addAnnotation(split[8].trim());
					}
				}
				if(lineCount % 10000 == 0)
					System.err.print(".");
				
			}
			System.err.println(args[gtf]+" had "+count);
			ss.close();
			ss = new SuperScanner(args[1]);
		}


		//ok now that we have loaded all of the GTF files lets compare them all against each other
		System.err.println();
		for(String loc: annotations.keySet())  //for each chr
		{
			System.err.print("Processing "+loc); 
			for(int gtf1: annotations.get(loc).keySet())
			{
				for(int gtf2: annotations.get(loc).keySet())
				{

					//lets go through and find the best score for each transcript
					for(Annotation a: annotations.get(loc).get(gtf1))
					{
						if(gtf1 != gtf2)
						{
							int bestScore = 0;
							for(Annotation b: annotations.get(loc).get(gtf2))
							{
								int score = a.overlaps(b);
								if(score > bestScore) bestScore = score;
							}
							scores.get(a.id)[gtf2]=bestScore;
						}
						else
							scores.get(a.id)[gtf2]=5;
					}
						
				}
			}
			System.err.println("\tDone");
		}
		//now for each unique transcript lets count up instances for all binary memberships in all gtfs

		String header = "Count";
		for(String arg: args)
		{
			header+="\t"+arg;
		}
		System.out.println(header);
		
		HashMap<String,Integer> counts = new HashMap<String,Integer>();
		for(Annotation a: byId.values())
		{
			int[] score = scores.get(a.id);
			String boolStr = "";
			String specific = ""+fileMap.get(a.id)+":";
			for(int i = 0; i < score.length; i++)
			{
				if(score[i] >0)
					boolStr+="1";
				else
					boolStr+="0";
			}
			specific+=boolStr;
			if(counts.get(boolStr) == null) counts.put(boolStr,1);
			else counts.put(boolStr,counts.get(boolStr)+1);
			if(counts.get(specific) == null) counts.put(specific,1);
			else counts.put(specific,counts.get(specific)+1);			
		}
		
		int maxValue = (int)Math.pow(2,args.length);
		int sum = 0;
		for(int i = 1; i < maxValue; i++)
		{
			String boolStr = Integer.toBinaryString(i);
			while(boolStr.length() < args.length)
				boolStr="0"+boolStr;
			if(counts.get(boolStr)!= null)
			{
				System.out.print(counts.get(boolStr));
				sum+=counts.get(boolStr);
			}
			else System.out.print("0");
			for(int j = 0; j < boolStr.length() ;j++)
			{

				if(boolStr.charAt(j)=='1')
					System.out.print("\tX");
				else System.out.print("\t");
			}
			System.out.println();

		}
		for(int gtf = 0; gtf < args.length; gtf++)
		{
			System.out.println("Specific counts for "+args[gtf]);
			for(int i = 1; i < maxValue; i++)
			{
				String boolStr = Integer.toBinaryString(i);
				while(boolStr.length() < args.length)
					boolStr="0"+boolStr;
				if(counts.get(gtf+":"+boolStr)!= null)
				{
					System.out.print(counts.get(gtf+":"+boolStr));
				}
				else System.out.print("0");
				for(int j = 0; j < boolStr.length() ;j++)
				{
					if(boolStr.charAt(j)=='1')
						System.out.print("\tX");
					else System.out.print("\t");
				}
				System.out.println();

			}
		}
		System.out.println("Sum:\t"+sum);
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new GTFVenn(args);
	}

}
