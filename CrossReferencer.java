import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


public class CrossReferencer {

	private class Peptide
	{
		String pep ="";
		HashMap<String,String> orfs = new HashMap<String,String>();  //orfTag-> nucseq //can only have one nuc seq per specific orf tag
		HashMap<String,String> orfStrands = new HashMap<String,String>(); //orfTag->true = positive strand
		HashMap<String,String> orfChr = new HashMap<String,String>(); //orfTag->true = positive strand
		HashMap<String,LinkedList<Integer>> orfExonCoords = new HashMap<String,LinkedList<Integer>>();
		public Peptide(String pep)
		{
			this.pep = pep;
		}
	}
	
	public CrossReferencer(String[] args) {
		
		LinkedList<Peptide> peptides = new LinkedList<Peptide>();  //all peptides detected by mass spec (arg0)
		HashSet<String> orfTags = new HashSet<String>() ; // all ORFs that had hits
		
		
		if(args.length < 4)
		{
			System.err.println("Usage: java CrossReferencer mass_spec_peptides.txt 3frame.fa nuc.fa orfs.gtf > peptides_with_data.txt");
			System.exit(1);
		}
		SuperScanner ss = new SuperScanner(args[0]);
		while(ss.hasMore())
		{
			String line = ss.getLine();
			if(line.length() > 0)
				peptides.add(new Peptide(line));
		}
		ss.close();
		ss = new SuperScanner(args[1]);
		System.err.println("Loaded "+peptides.size()+" peptides from "+args[0]);
		String tag = null;
		String seq = "";
		System.err.println("Searching 3frame for peptides");
		while(ss.hasMore())
		{
			String line = ss.getLine();
			if(line.startsWith(">"))
			{
				if(tag != null)
				{
					//save the tag if it contains one of the peptides
					for(Peptide pep: peptides)
					{
						if(seq.contains(pep.pep))
						{
							System.err.println("ORF: "+tag+" contained: "+pep.pep);
							pep.orfs.put(tag, null); // put an empty nuc ref to be populated later
							orfTags.add(tag);
						}
					}

				}
				tag = line.substring(1);
				seq = "";
			}
			else
				seq+= line;
		}
		System.err.println("Done searching the peptide dataset...");
		if(tag != null)
		{
			//save the tag if it contains one of the peptides
			for(Peptide pep: peptides)
			{
				if(seq.contains(pep.pep))
				{
					pep.orfs.put(tag, null); // put an empty nuc ref to be populated later
					orfTags.add(tag);
				}
			}

		}
		ss.close();
		ss = new SuperScanner(args[2]);
		System.err.println("Searching the nuc file for corresponding nuc sequence");
		while(ss.hasMore())
		{
			String line = ss.getLine();
			if(line.startsWith(">"))
			{
				if(tag != null)
				{
					//save the tag if it contains one of the peptides
//					if(tag.contains("CUFF.31297.1+chr4:62658443-62658502_F:2_P:6"))
//						System.err.println();
					for(Peptide pep: peptides)
					{
						HashMap<String,String> temp = new HashMap<String,String>();
						for(String orf: pep.orfs.keySet())
						{
							//this is the name of the orf of interest... lets check to see if it is contained in the nuc
							if(tag.contains(orf))
							{
								temp.put(orf, seq); 
								System.err.println("Found nuc entry for :"+orf);
							}
							else
								temp.put(orf, pep.orfs.get(orf));
						}
						pep.orfs = temp; //update at the end to avoid concurrent modification of the hashmap						
					}
				}
				tag = line.substring(1);
				seq = "";
			}
			else
				seq+= line;
		}
		if(tag != null)
		{
			//save the tag if it contains one of the peptides
			for(Peptide pep: peptides)
			{
				HashMap<String,String> temp = new HashMap<String,String>();
				for(String orf: pep.orfs.keySet())
				{
					//this is the name of the orf of interest... lets check to see if it is contained in the nuc
					if(tag.contains(orf))
					{
						temp.put(orf, seq);
						System.err.println("Found nuc entry for "+orf);
					}
					else
						temp.put(orf, pep.orfs.get(orf));
				}
				pep.orfs = temp; //update at the end to avoid concurrent modification of the hashmap						
			}
		}
		ss.close();
		//finally lets go through the GTF file searching for relevant coords
		System.err.println("Grabbing relevant coordinates from the GTF file and outputing to STDOUT");
		System.err.println("The output will be: ORFcoord[tab]ORFtag[tab]peptide[tab]transcriptSeq[tab]exonCoords");
		ss = new SuperScanner(args[3]);
		while(ss.hasMore())
		{
			String[] split = ss.getLine().split("\t");
			String id = split[8].substring(split[8].indexOf("transcript_id")+15,split[8].indexOf('"', split[8].indexOf("transcript_id")+15));
			if(orfTags.contains(id))
			{
				//this entry is of potential use!
				if(split[2].compareTo("exon")==0)
				{
					//lets save the strandedness and exon info
					for(Peptide pep: peptides)
					{
						for(String peptag: pep.orfs.keySet())
						{
							if(peptag.compareTo(id)==0)
							{
								System.err.println("Found gtf entry for "+peptag);
								//we found the orf of interest
								pep.orfStrands.put(peptag, split[6]);
								pep.orfChr.put(peptag,split[0]);
								if(pep.orfExonCoords.get(peptag)== null) pep.orfExonCoords.put(peptag, new LinkedList<Integer>());
								pep.orfExonCoords.get(peptag).add(Integer.parseInt(split[3]));
								pep.orfExonCoords.get(peptag).add(Integer.parseInt(split[4]));
							}
						}
					}
				}
					
			}
		}
		for(Peptide pep: peptides)
		{
			for(String orf: pep.orfs.keySet())
			{
				//lets output the hit!
				System.out.print(pep.orfChr.get(orf)+":"+pep.orfExonCoords.get(orf).getFirst()+"-"+pep.orfExonCoords.get(orf).getLast()+" strand="+pep.orfStrands.get(orf)+"\t"+orf+"\t"+pep.pep+"\t"+pep.orfs.get(orf)+"\t");
				for(int i = 0; i < pep.orfExonCoords.get(orf).size(); i+=2)
					System.out.print(pep.orfExonCoords.get(orf).get(i)+","+pep.orfExonCoords.get(orf).get(i+1)+";");
				System.out.println();
			}
		}
		System.err.println("Finished cross referencing the databases!");
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new CrossReferencer(args);
	}

}
