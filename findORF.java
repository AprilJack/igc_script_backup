import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;



public class findORF {

	public findORF(String[] args) {
		HashMap<String,char[]> genomeMap = new HashMap<String,char[]>();
		HashMap<String,String> nucToPep = buildNucToPep();
		try {
			Scanner s = new Scanner(new File(args[0]));
			StringBuilder sb = new StringBuilder();
			String reading = null;
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				if(line.startsWith(">"))
				{
					if(reading != null)
					{
						char[] chr = new char[sb.length()];
						sb.getChars(0, sb.length(), chr, 0);
						genomeMap.put(reading,chr);
					}
					reading = line.substring(1);
					sb = new StringBuilder();
				}
				else
				{
					sb.append(line);
				}
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for(int i = 1; i < args.length; i++)
		{
			try{
				String chr = args[i].split(":")[0];
				int pos = Integer.parseInt(args[i].split(":")[1]);
				boolean negative = (args[i].split(":")[2].compareTo("-")==0);
				int p = pos;
				String result = "";
				char[] seq = genomeMap.get(chr);
				if(!negative)
				{
					while(p < seq.length)
					{
						String trip = ""+seq[p]+seq[p+1]+seq[p+2];
						String pep = nucToPep.get(trip);
						if(pep.compareTo("X")!= 0)
							result+=pep;
						else break;
					}
				}
			}catch(Exception e){System.err.println("Failed to process "+args[i]);}
		}
	}
	
	private HashMap<String, String> buildNucToPep() {
		// TODO Auto-generated method stub
		HashMap<String,String> map = new HashMap<String,String>();
		map.put("TTT","F");map.put("TCT","S");map.put("TAT","Y");map.put("TGT","C");
		map.put("TTC","F");map.put("TCC","S");map.put("TAC","Y");map.put("TGC","C");
		map.put("TTA","L");map.put("TCA","S");map.put("TAA","X");map.put("TGA","X");
		map.put("TTG","L");map.put("TCG","S");map.put("TAG","X");map.put("TGG","W");
		map.put("CTT","L");map.put("CCT","P");map.put("CAT","H");map.put("CGT","R");
		map.put("CTC","L");map.put("CCC","P");map.put("CAC","H");map.put("CGC","R");
		map.put("CTA","L");map.put("CCA","P");map.put("CAA","Q");map.put("CGA","R");
		map.put("CTG","L");map.put("CCG","P");map.put("CAG","Q");map.put("CGG","R");
		map.put("ATT","I");map.put("ACT","T");map.put("AAT","N");map.put("AGT","S");
		map.put("ATC","I");map.put("ACC","T");map.put("AAC","N");map.put("AGC","S");
		map.put("ATA","I");map.put("ACA","T");map.put("AAA","K");map.put("AGA","R");
		map.put("ATG","M");map.put("ACG","T");map.put("AAG","K");map.put("AGG","R");
		map.put("GTT","V");map.put("GCT","A");map.put("GAT","D");map.put("GGT","G");
		map.put("GTC","V");map.put("GCC","A");map.put("GAC","D");map.put("GGC","G");
		map.put("GTA","V");map.put("GCA","A");map.put("GAA","E");map.put("GGA","G");
		map.put("GTG","V");map.put("GCG","A");map.put("GAG","E");map.put("GGG","G");
		HashMap<String,String> newMap = new HashMap<String,String>();
		for(String key: map.keySet())
		{
			String copy = new String(key);
			String TtoU=copy.replace('T', 'U');
			newMap.put(TtoU,map.get(key));
			newMap.put(key,map.get(key));
		}
		return newMap;	
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 2)
		{
			System.err.println("Usage: findORF genome.fa chr:pos:- [chr:pos:+] ...");
			System.err.println("Searches upstream and downstream of the position specified until it finds a start and stop codon, respectively (assumes in frame). Specify - for reverseComplement search. Outputs in fasta format.");
			System.exit(1);
		}
		new findORF(args);
	}

}
