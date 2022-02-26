import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;


public class SNPcompare {

	HashMap<String,Annotation> byId = new HashMap<String,Annotation>();
	HashMap<String,LinkedList<Annotation>> byChr = new HashMap<String,LinkedList<Annotation>>();
	HashMap<String,LinkedList<String>> SNPs = new HashMap<String,LinkedList<String>>();
	
	public SNPcompare(String[] args) {
		if(args.length < 3)
		{
			System.err.println("Usage: SNPcompare Genes.gtf File1.allSNPs.txt File2.allSNPs.txt");
			System.exit(1);
		}
		SuperScanner ss = new SuperScanner(args[0]);
		System.err.print("Loading "+args[0]);
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
				if(ann.compareTo("exon")==0 )
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
		ss = new SuperScanner(args[1]);
		System.err.print("Loading "+args[1]);
		lineCount = 0;
		while(ss.hasMore())
		{
			//load in one line at a time
			try{
				String line = ss.getLine();
				String[] split = line.split("\t");
				String chrStr = split[1];
				int pos =  Integer.parseInt(split[2])+1;
				if(SNPs.get(chrStr+":"+pos)==null)SNPs.put(chrStr+":"+pos, new LinkedList<String>());
				SNPs.get(chrStr+":"+pos).add(split[5]);
			}catch(Exception e)
			{
				
			}
			if(lineCount % 10000 == 0)
				System.err.print(".");
			lineCount++;
		}
		System.err.println();
		ss = new SuperScanner(args[2]);
		System.err.print("Seeing how "+args[2]+" is different and outputing gene summary to stdout");
		lineCount = 0;
		while(ss.hasMore())
		{
			try{
				String line = ss.getLine();
				String[] split = line.split("\t");
				String chrStr = split[1];
				int pos =  Integer.parseInt(split[2])+1;
				//lets check if the SNP was already there
				if(SNPs.get(chrStr+":"+pos)==null)SNPs.put(chrStr+":"+pos, new LinkedList<String>());
				SNPs.get(chrStr+":"+pos).add("!"+split[5]);//second one has a !

			}
			catch(Exception e)
			{
				
			}
			if(lineCount % 10000 == 0)
				System.err.print(".");
			lineCount++;
		}
		System.err.println("Outputing the results in a gene-wise context");
		{
			HashMap<String,Integer> indecies = new HashMap<String,Integer>();
			indecies.put("A/C",0);
			indecies.put("A/G",1);
			indecies.put("A/T",2);
			indecies.put("C/A",3);
			indecies.put("C/G",4);
			indecies.put("C/T",5);
			indecies.put("G/A",6);
			indecies.put("G/C",7);
			indecies.put("G/T",8);
			indecies.put("T/A",9);
			indecies.put("T/C",10);
			indecies.put("T/G",11);
			String header = "Gene\tChr\tPLUS\tLength(bp)\tSNP locations(bps)\tsame\tdifferent\tnovel\tmissing\ttotal";
			for(String mut: indecies.keySet())
				header+="\t1:"+mut;

			for(String mut: indecies.keySet())
				header+="\t2:"+mut;
			System.out.println(header);
			for(String chrStr: byChr.keySet())
			{
				boolean plus = chrStr.startsWith("+");
				String chr = chrStr.substring(1,chrStr.length());
				if(byChr.get(chrStr) != null)
				{

					for(Annotation a: byChr.get(chrStr))
					{
						int novel = 0;
						int missing = 0;
						int same = 0;
						int different = 0;
						int total = 0;
						int[] mold = new int[12];   
						int[] mnew = new int[12];
						int length = 0;
						String positions ="";
						//lets scan the entire annotation
						for(Region r: a.getRegions())
						{
							length+=1+r.end-r.start;
							for(int i = r.start; i<=r.end; i++)
							{
								if(SNPs.get(chr+":"+i) != null)
								{
									//lets figure out if it is the same or different
									if(SNPs.get(chr+":"+i).size()==1)
									{
										//either a novel or missing SNP
										if(SNPs.get(chr+":"+i).get(0).startsWith("!"))
										{
											novel++;
											mnew[indecies.get(strandedAnnotation(SNPs.get(chr+":"+i).get(0).substring(1),plus))]++;
											positions+=" "+i+"(novel)";
										}
										else {
											missing++;
											mold[indecies.get(strandedAnnotation(SNPs.get(chr+":"+i).get(0),plus))]++;
											positions+=" "+i+"(missing)";
										}
									}
									else if (SNPs.get(chr+":"+i).size()==2)
									{
										if(SNPs.get(chr+":"+i).get(1).substring(1).compareTo(SNPs.get(chr+":"+i).get(0))==0)
										{
											same++;
											positions+=" "+i+"(same)";
										}
										else
										{
											different++;
											positions+=" "+i+"(diff)";
										}
										mold[indecies.get(strandedAnnotation(SNPs.get(chr+":"+i).get(0),plus))]++;
										mnew[indecies.get(strandedAnnotation(SNPs.get(chr+":"+i).get(1).substring(1),plus))]++;
									}
									
									total++;
								}
							}
						}
						String out = a.id+"\t"+chr+"\t"+plus+"\t"+length+"\t"+positions+"\t"+same+"\t"+different+"\t"+novel+"\t"+missing+"\t"+total;
						for(String mut: indecies.keySet())
							out+="\t"+mold[indecies.get(mut)];
						for(String mut: indecies.keySet())
							out+="\t"+mnew[indecies.get(mut)];
						System.out.println(out);
						
					}
				}
			}
		}
	}
	
	private String strandedAnnotation(String mut, boolean plus)
	{
		if(plus) return mut;
		else
		{
			char from = mut.split("/")[0].charAt(0);
			char to = mut.split("/")[1].charAt(0);
			if(from == 'A') from = 'T';
			else if(from =='T') from = 'A';
			else if(from =='C') from = 'G';
			else if(from =='G') from = 'C';
			if(to == 'A') to = 'T';
			else if(to =='T') to = 'A';
			else if(to =='C') to = 'G';
			else if(to =='G') to = 'C';
			return from+"/"+to;
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new SNPcompare(args);
	}

}
