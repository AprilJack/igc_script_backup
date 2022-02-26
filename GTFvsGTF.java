import java.util.HashMap;
import java.util.LinkedList;


public class GTFvsGTF {

	HashMap<String,LinkedList<Annotation>> gtf1 = new HashMap<String,LinkedList<Annotation>>();
	HashMap<String,LinkedList<Annotation>> gtf2 = new HashMap<String,LinkedList<Annotation>>();
	HashMap<Annotation,Annotation> closest = new HashMap<Annotation,Annotation>();
	HashMap<String, Annotation> byId1 = new HashMap<String,Annotation>();
	HashMap<String, Annotation> byId2 = new HashMap<String,Annotation>();
	HashMap<String, LinkedList<String>> gtfLines = new HashMap<String,LinkedList<String>>();
	int totalGTF1 = 0;
	int totalGTF2 = 0;
	String[] codes = {"UN","SO","SM","AO","IEqual","Equal"};
	boolean debug = false;
	
	public GTFvsGTF(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: GTFvsGTF gtf1.txt gtf2.txt [size-cutoff]  > gtf1_annotated.gtf\n"
					+ "Compares two GTF files and annotates gtf1 with gtf2. Counts total transcripts, overlapping, \n"
					+ "and partially overlapping UTRs, and partially overlapping missing-exons.\n"
					+ "If a GTF has CDSs then only coding transcripts will be used. Annotated.gtf is gtf1 plus a match code.\n"
					);
			System.exit(1);
		}
		SuperScanner ss = new SuperScanner(args[0]);
		boolean GTF1HasCDS = false;
		boolean GTF2HasCDS = false;
		System.err.println("Loading "+args[0]);
		int lineCount = 0;
		while(ss.hasMore())
		{
			lineCount++;
			String line = ss.getLine();
			String[] split = line.split("\t");
			if(split.length > 8 && split[2].length() >0)
			{
				String chrStr = split[0]+split[6];
				String ann = split[2];
				if(ann.compareTo("exon")==0)
				{
					int start = Integer.parseInt(split[3].trim());
					int end = Integer.parseInt(split[4].trim());
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					String key = chrStr;
					if(gtf1.get(key) == null) gtf1.put(key, new LinkedList<Annotation>());
					if(byId1.get(id) == null)
					{
						byId1.put(id,new Annotation(chrStr,id));
						gtf1.get(key).add(byId1.get(id));
						gtfLines.put(id, new LinkedList<String>());
						totalGTF1++;
					}
					byId1.get(id).addRegion(start, end);
					byId1.get(id).addAnnotation(split[8].trim());
					gtfLines.get(id).add(line);
				}
				else if(ann.compareTo("CDS")==0)
				{
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					String key = chrStr;
					if(gtf1.get(key) == null) gtf1.put(key, new LinkedList<Annotation>());  //if the chr doesn't yet exist
					if(byId1.get(id) == null)  //if the id doesn't yet exist
					{
						byId1.put(id,new Annotation(chrStr,id));    //add the new annotation
						gtf1.get(key).add(byId1.get(id));   //add the annotation to the gtf1 list
						byId1.get(id).coding=true;  //set it to coding
						gtfLines.put(id, new LinkedList<String>());
						totalGTF1++;
						GTF1HasCDS=true;  
					}
					else  //already exists	
					{
						byId1.get(id).coding = true;
						GTF1HasCDS =true;
					}
					gtfLines.get(id).add(line);
				}
			}
			if(lineCount % 10000 == 0)
				System.err.print(".");
			
		}
		ss.close();
		ss = new SuperScanner(args[1]);
		
		System.err.println();
		System.err.println("Loading "+args[1]);
		while(ss.hasMore())
		{
			lineCount++;
			String line = ss.getLine();
			if(line == null) break;
			String[] split = line.split("\t");
			if(split.length > 8 && split[2].length() >0)
			{
				String chrStr = split[0]+split[6];
				String ann = split[2];
				if(ann.compareTo("exon")==0)
				{
					int start = Integer.parseInt(split[3].trim());
					int end = Integer.parseInt(split[4].trim());
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					String key = chrStr;
					if(gtf2.get(key) == null) gtf2.put(key, new LinkedList<Annotation>());
					if(byId2.get(id) == null)
					{
						byId2.put(id,new Annotation(chrStr,id));
						gtf2.get(key).add(byId2.get(id));

						totalGTF2++;
					}
					byId2.get(id).addRegion(start, end);
					byId2.get(id).addAnnotation(split[8].trim());

				}
				else if(ann.compareTo("CDS")==0)
				{
					String id = split[8].replaceAll("\"", "");
					id=id.replaceAll("\"", "");
					id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
					String key = chrStr;
					if(gtf2.get(key) == null) gtf2.put(key, new LinkedList<Annotation>());
					if(byId2.get(id) == null)
					{
						byId2.put(id,new Annotation(chrStr,id));
						gtf2.get(key).add(byId2.get(id));
						byId2.get(id).coding=true;
						totalGTF2++;
						GTF2HasCDS=true;
					}
					else
					{
						GTF2HasCDS=true;
						byId2.get(id).coding =true;
					}

				}
			}
			if(lineCount % 10000 == 0)
				System.err.print(".");
		}
		ss.close();
		int size = Integer.MAX_VALUE;
		if(args.length > 2)
			size = Integer.parseInt(args[2]);
		//ok now that we have loaded all of the GTF files lets compare them all against each other
		System.err.println();
		int[] scoreBins = new int[6];
		int totalScores = 0;
		for(String loc: gtf1.keySet())  //for each chr
		{
			System.err.println("Processing "+loc); 
			int count = 0;
			for(Annotation a: gtf1.get(loc))  //for each gtf1 on the chr
			{
				int bestScore = 0;
				Annotation bestB = null;
				if(gtf2.get(loc) != null)
				{
					for(Annotation b: gtf2.get(loc))  //for each annotation on the same chr of gtf2
					{
						int score = a.overlaps(b);
						if(score > bestScore)
						{
							bestScore = score;
							bestB = b;
						}
					}
				}
				scoreBins[bestScore]++;
				totalScores++;
				count++;
				if(!debug)
				{
					for(String line: gtfLines.get(a.id))
					{
						if(bestB != null)
							System.out.println(line+" Code "+codes[bestScore]+"; Closest: "+bestB.id);
						else
							System.out.println(line+" Code "+codes[bestScore]);
					}
					
				}
				if(count% 100 == 0)
					System.err.print(".");
			}
		}
		System.err.println();
		/* 1 = some exons are overlapped
				 * 2 = some exons are perfectly overlapped
				 * 3 = all exons are overlapped
				 * 4 = internal exons are perfectly matched
				 * 5 = all exons are perfectly matched (identical) */
		System.err.println("Reference\tVS\tSum Ref\tSum Vs\tTotal Scores\tNo Match\tsome overlaped\tsome matched\tall overlapped\tinternal matched\tperfect match");
		System.err.println(args[0]+"\t"+args[1]+"\t"+totalGTF1+"\t"+totalGTF2+"\t"+totalScores+"\t"+scoreBins[0]+"\t"+scoreBins[1]+"\t"+scoreBins[2]+"\t"+scoreBins[3]+"\t"+scoreBins[4]+"\t"+scoreBins[5]);
		System.err.println("Done!");
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new GTFvsGTF(args);
	}

}
