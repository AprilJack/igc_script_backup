import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;


public class PeakExtractor {

	static HashMap<String,byte[]> chrs = new HashMap<String,byte[]>();
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 3)
		{
			System.err.println("Usage: peakExtractor hg19.fa Riboseq.bedGraph peak_threhsold");
			System.exit(1);;
		}
		double threshold = Double.parseDouble(args[2]);
		
		try {
			Scanner s = new Scanner(new File(args[0]));
			PrintWriter fasta = new PrintWriter(args[1].substring(0, args[1].lastIndexOf('.'))+".fasta");  //fasta 
			PrintWriter threeframe = new PrintWriter(args[1].substring(0, args[1].lastIndexOf('.'))+".3frame");  //fasta 3frame
			PrintWriter bed = new PrintWriter(args[1].substring(0, args[1].lastIndexOf('.'))+".bed"); //bed of alignment
			StringBuilder sb = null;
			String chr = null;
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				if(line.startsWith(">"))
				{
					if(chr != null)
					{
						chrs.put(chr, sb.toString().getBytes());
						System.err.println("Loaded chr: "+chr+" with "+chrs.get(chr).length+" bases. Memory free: "+(Runtime.getRuntime().freeMemory()/(1024*1024))+" MB / "+(Runtime.getRuntime().totalMemory()/(1024*1024))+" MB");
					}
					chr = line.substring(1);
					sb = new StringBuilder(100000000);
				}
				else
				{
					sb.append(line);
				}
			}
			s.close();
			s = new Scanner(new File(args[1]));
			String seq = null;
			double intensity =0;
			String chrName =null;
			int seqStart = 0;
			int count = 0;
			
			int seg_count = 0;
			boolean forward = true;
			boolean in_peak = false;
			long peakLengths = 0;
			long segLengths = 0;
			String lastEnd = null;
			System.err.println();
			bed.println("track name="+args[1].substring(args[1].lastIndexOf('/'))+" description=\"bedGraph peaks with threshold "+args[2]+"\" useScore=1");
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				if(line.startsWith("track type="))
				{
					if(seq != null)
					{
						if(in_peak)
						{
							/*lets save the peak 3frame translated
							for(int i = 0; i < 3; i++)
							{
								String[] translated = getFrameTranslated(seq, i).split("X");
								int seqThusFar = 0;
								for(int j = 0; j < translated.length; j++)
								{
									if(translated[j].length() > 9 && translated[j].indexOf('?')==-1)
									{
										fasta.println(">Peak "+count+" FWD: "+forward+" frame "+(i+1)+" segment "+j);						
										fasta.println(translated[j]);
										String strand = "+"; if(!forward) strand = "-";
										bed.println(chrName+"\t"+(seqStart+i+seqThusFar+1)+"\t"+(seqStart+i+seqThusFar+translated[j].length()*3+1)+"\t"+
												"Peak"+count+"_frame"+(i+1)+"_seg"+j+"\t"+(int)(Math.min(1000,intensity/seq.length()))+"\t"+strand);
										seg_count++;
										segLengths+=translated[j].length();
									}
									if(j ==0)
										seqThusFar+=(translated[j].length()*3);
									else seqThusFar+=(3+translated[j].length()*3);
								}
							}
							*/
							in_peak = false;
						}
						forward =false;
					}
				}
				else
				{
					
					String[] split = line.split("\\s+");
					double value = Double.parseDouble(split[3]);
					chrName = split[0];
					if(value > threshold && !in_peak )
					{
						if(forward)
							seq = getSequence(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]));
						else
							seq = revComplement(getSequence(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2])));
						in_peak = true;
						seqStart = Integer.parseInt(split[1]);
						intensity +=value*(Integer.parseInt(split[2])-Integer.parseInt(split[1]));
						lastEnd = split[2];
					}
					else if(value > threshold && in_peak &&(lastEnd == null || lastEnd.compareTo(split[1])==0))
					{
						if(forward)
							seq += getSequence(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]));
						else
							seq = revComplement(getSequence(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2])))+seq;
						lastEnd = split[2];
					}
					else if(in_peak &&((value < threshold)|| (lastEnd != null &&lastEnd.compareTo(split[1])!=0))) //fell below or missing values
					{
						if(seq.length() > 30)  //length check
						{
							peakLengths+=(seq.length()-seq.length()%3);
							count++;
							if(count%1000 ==0) System.err.print(".");
							//we are done with this seq!
							fasta.println(">Peak "+count+" FWD: "+forward);							
							fasta.println(seq);
							//lets save the peak 3frame translated and extended up and downstream
							for(int i = 0; i < 3; i++)
							{
								String translatedStr = "";
								translatedStr=getFrameTranslated(seq, i);
								if(translatedStr.lastIndexOf('X')==translatedStr.indexOf('X')&&!translatedStr.contains("?"))
								{
									if(translatedStr.indexOf('X')!=-1)
										translatedStr=translatedStr.substring(0,translatedStr.indexOf('X'));
									
									if(translatedStr.length() > 9 && chrs.get(chrName).length > (seqStart+i+translatedStr.length()*3))
									{
										//String fullORF = scanUpAndDown(chrName,seqStart,seqStart+seq.length(),i,forward);
										String strand = "+"; if(!forward) strand = "-";
										
										if(forward)
										{
											int frame = (i+seqStart)%3;
											threeframe.println(">Peak"+count+"_frame"+(frame+1)+"_"+strand);						
											threeframe.println(translatedStr);
											bed.println(chrName+"\t"+(seqStart+i)+"\t"+(seqStart+i+translatedStr.length()*3)+"\t"+
												"Peak"+count+"_frame"+(frame+1)+"\t"+(int)(Math.min(1000,500+intensity/seq.length()))+"\t"+strand);
										}
										else
										{
											int leftOffset = seq.length()-translatedStr.length()*3-i;
											int frame = (i+seqStart+1)%3;
											threeframe.println(">Peak"+count+"_frame"+(frame+1)+"_"+strand);						
											threeframe.println(translatedStr);
											bed.println(chrName+"\t"+(seqStart+leftOffset)+"\t"+(seqStart+leftOffset+translatedStr.length()*3)+"\t"+
													"Peak"+count+"_frame"+(frame+1)+"\t"+(int)(Math.min(1000,500+intensity/seq.length()))+"\t"+strand);
										}
										seg_count++;
										segLengths+=translatedStr.length()*3;
									}
								}
							}
						}
						in_peak = false;
					}
				}
			}
			s.close();
			fasta.close();
			threeframe.close();
			bed.close();
			System.err.println(String.format("Total peaks:%d, Avg. Peak Length: %3.1f, saved segments (across 3 frames): %d, avg. seg length: %3.1f", count,((double)(peakLengths/count)),seg_count,((double)(segLengths/seg_count))));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static String scanUpAndDown(String chr, int start, int end, int offset, boolean forward)
	{
		String seq = "";
		if(forward)
		{
			int pos = start+offset;
			while(true)
			{
				try{
					String AA = seqToAA(getSequence(chr,pos-3,pos));
					if(AA.compareTo("M")==0)
					{
						seq="M"+seq;
						break;
					}
					else if(AA.compareTo("X")==0)
					{
						break;
					}
					else
						seq=AA+seq;
					pos-=3;
				}
				catch(Exception e){seq="!"+seq; break;}
			}
			pos = start+offset;
			while(true)
			{
				//scan down
				try{
					String AA = seqToAA(getSequence(chr,pos,pos+3));
					if(AA.compareTo("X")==0)
					{
						break;
					}
					else
						seq=seq+AA;
					pos+=3;
				}
				catch(Exception e){seq=seq+"!"; break;}
			}
		}
		else //reverse
		{
			int pos = end-offset;
			while(true)
			{
				try{
					String AA = seqToAA(revComplement(getSequence(chr,pos,pos+3)));
					if(AA.compareTo("M")==0)
					{
						seq="M"+seq;
						break;
					}
					else if(AA.compareTo("X")==0)
					{
						break;
					}
					else
						seq=AA+seq;
					pos+=3;
				}
				catch(Exception e){seq="!"+seq; break;}
			}
			pos = end-offset;
			while(true)
			{
				//scan down
				try{
					String AA = seqToAA(revComplement(getSequence(chr,pos-3,pos)));
					if(AA.compareTo("X")==0)
					{
						break;
					}
					else
						seq=seq+AA;
					pos-=3;
				}
				catch(Exception e){seq=seq+"!"; break;}
			}
		}
		return seq;
		
	}
	
	private static String getSequence(String chr, int start,
			int stop) {
		char[] chars = new char[(int) (stop-start)];
		if(chrs.get(chr) != null)
		{
			byte[] bytes = chrs.get(chr);
			try{
				for(int i = start; i < stop; i++)
					chars[i-start]=(char) bytes[i];
			}catch(Exception e){}
		}
		return new String(chars).toUpperCase();
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
