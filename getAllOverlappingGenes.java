import java.util.HashMap;
import java.util.LinkedList;


public class getAllOverlappingGenes {

	
	HashMap<String,HashMap<String,Integer>> transcripts = new HashMap<String,HashMap<String,Integer>>();
	
	
	public getAllOverlappingGenes(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Outputs a list of all genes overlapping the regions of interest (assume TSS as overlap point and both strands)\n"
					+ "\nUsage: getAllOverlappingGenes peak_files.txt annotation.gtf [distance]> geneList.txt");
			System.exit(1);
		}
		int distance = 0;
		if (args.length > 2)
			distance = Integer.parseInt(args[2]);
		System.err.println("Processing "+args[1]+" using "+args[0]);
		SuperScanner gtfScanner = new SuperScanner(args[1]);
		while(gtfScanner.hasMore())
		{
			String line = gtfScanner.getLine();
			String[] split = line.split("\t");
			if(split.length < 3) continue;
			String chr = split[0];
			String ann = split[2];
			if(ann.compareTo("exon")==0)
			{
				int start = Integer.parseInt(split[3].trim());
				int end = Integer.parseInt(split[4].trim());
				String id = split[8].replaceAll("\"", "");
				id=id.replaceAll("\"", "");
				id=id.substring(id.indexOf("transcript_id")+14,id.indexOf(";", id.indexOf("transcript_id")+14));
				if(transcripts.get(chr) == null) transcripts.put(chr,new HashMap<String,Integer>());
				String transcript = new String(line.substring(line.indexOf("transcript_id")+15));
				transcript = transcript.substring(0,transcript.indexOf("\""));
				if(transcripts.get(chr).get(transcript) == null)
				{
					if(split[6].compareTo("+")==0)
						transcripts.get(chr).put(transcript,start);
					else
						transcripts.get(chr).put(transcript,end);
				}
				else
				{
					if(split[6].compareTo("+")==0)
						transcripts.get(chr).put(transcript,Math.min(transcripts.get(chr).get(transcript),start));
					else
						transcripts.get(chr).put(transcript,Math.max(transcripts.get(chr).get(transcript),end));
				}
			}
		}
		gtfScanner.close();
		SuperScanner posScanner = new SuperScanner(args[0]);
		int totalRegions = 0;
		int totalGenes = 0;
		while(posScanner.hasMore())
		{
			String line = posScanner.getLine();
			String[] split = line.split("\t");
			try{
				String chr = split[1];
				int min = Integer.parseInt(split[2]);
				int max = Integer.parseInt(split[3]);
				totalRegions++;
				if(transcripts.get(chr) != null)
				{
					for(String gene: transcripts.get(chr).keySet())
					{
						int tss = transcripts.get(chr).get(gene);
						if(tss+distance >= min && tss-distance <= max)
						{
							System.out.println(gene+"\t"+split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]);
							totalGenes++;
						}
					}
				}
			}
			catch(Exception e){}
		}
		System.err.println("Finished processing "+totalRegions+" regions which had "+totalGenes+" genes.");
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new getAllOverlappingGenes(args);
	}

}
