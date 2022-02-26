import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.omg.CORBA.INITIALIZE;


public class positionsToSeq {

	HashMap<String,byte[]> chrs = new HashMap<String,byte[]>();
	
	public positionsToSeq(String positionsFile, String sequencesFile, int chrCol, String delim, int scoreThreshold) {
	
			FastScanner s = new FastScanner(sequencesFile);
			StringBuilder sb = new StringBuilder(300000000);
			String chr = null;
			String prefix = null;
			nucToPep =buildNucToPep();
			while(s.open)
			{
				String line = s.getLine();
				if(line.startsWith(">"))
				{
					if(chr != null)
					{
						
						chrs.put(chr,sb.toString().getBytes());
						sb = new StringBuilder(300000000);
						System.err.println("Loaded: "+chr+" with "+chrs.get(chr).length);
						chr = line.substring(1);
						
						prefix=getCommonPrefix(prefix,chr);

					}
					else
					{
						chr = line.substring(1);
						prefix = new String(chr);
						
					}
				}
				else
				{
					sb.append(line);
				}
			}
			chrs.put(chr,sb.toString().getBytes());
			HashMap<String,byte[]> temp = new HashMap<String,byte[]>();
			for(String key: chrs.keySet())
			{
				temp.put(key.substring(prefix.length()),chrs.get(key));
				temp.put("\""+key.substring(prefix.length())+"\"",chrs.get(key));
			}
			temp.putAll(chrs);
			chrs = temp;
			sb = null;
			System.err.println("Done reading: "+sequencesFile);
			s = new FastScanner(positionsFile);
			{
				//lets assume we have three consecutive columns starting the chrName
				int count = 0;
				while(s.max >0)
				{
					if(count++ %10000 == 0) System.err.print("*");
					String line = s.getLine();
					String[] split = line.split(delim);			
					try{
						if(chrs.get(split[chrCol])!= null)
						{

							int start = Integer.parseInt(split[chrCol+1]);
							int end = Integer.parseInt(split[chrCol+2]);
							int strand = Integer.parseInt(split[chrCol+3]);
							String sequence = getSeq(chrs.get(split[chrCol]),start-1,end-1,strand);
							
							if(!sequence.contains("?") && Integer.parseInt(split[chrCol+4]) >= scoreThreshold )
							{
								if(strand ==1)
									System.out.println(">S+"+split[chrCol]+":"+split[chrCol+1]+"-"+split[chrCol+2]);
								else
									System.out.println(">S-"+split[chrCol]+":"+split[chrCol+1]+"-"+split[chrCol+2]);
								System.out.println(sequence.substring(0,sequence.length()-1));
							}
						}
					}catch(Exception e){}
				}
				System.err.println();
				System.err.println("Done!");
			}
		
	}
	private String getCommonPrefix(String s1, String s2)
	{
		int i = 0;
		while(i < s1.length() && i < s2.length())
		{
			if(s1.charAt(i)!= s2.charAt(i))
				return s1.substring(0,i);
			i++;
		}
		if(s1.length() < s2.length())
			return s1;
		else return s2;

	}
	
	private String getSeq(byte[] buff, int start, int end, int strand)
	{
		if(end < start) 
		{
			int temp = end;
			end = start;
			start = temp;
		}
		byte[] buf2 = new byte[end-start+1];

		for(int i = start; i <= end; i++)
		{
			buf2[i-start] = buff[i];
		}
		String result = new String(buf2);
		if(strand == 1)
			return getFrameTranslated(result, 0).toUpperCase();
		else return getFrameTranslated(revComplement(result), 0).toUpperCase();
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new positionsToSeq(args[0], args[1],Integer.parseInt(args[2])-1,args[3],Integer.parseInt(args[4]));
	}
	
	private static HashMap<String, String> buildNucToPep() {
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
			newMap.put(key.toLowerCase(),map.get(key).toLowerCase());
			newMap.put(TtoU.toLowerCase(), map.get(key).toLowerCase());
		}
		return newMap;	
	}
	
	private static String getFrameTranslated(String seq, int frame)
	{
		String shift = seq.substring(frame);
		String trans = seqToAA(shift);
		return trans;
	}
	
	private static String revComplement(String seq)
	{
		String res = "";
		for(int i = seq.length()-1; i >= 0; i--)
		{
			if(seq.charAt(i)=='A') res+="T";
			else if(seq.charAt(i)=='T') res+="A";
			else if(seq.charAt(i)=='G') res+="C";
			else if(seq.charAt(i)=='C') res+="G";
			else res+=seq.charAt(i); //for now
		}
		return res;
	}
	static HashMap<String,String> nucToPep = buildNucToPep(); 
	
	private static String seqToAA(String seq)
	{
		String res = "";
		char[] arr = seq.toCharArray();
		for(int i = 0; i <= seq.length()-3; i+=3 )
		{
			String trip = ""+arr[i]+arr[i+1]+arr[i+2];
			String aa = "?";
			if(nucToPep.get(trip) != null)
				aa = nucToPep.get(trip);
			res+=aa;
		}
		return res;
	}

}
