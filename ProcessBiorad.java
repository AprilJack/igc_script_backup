import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


public class ProcessBiorad {

	HashMap<String,String[]> reads = new HashMap<String,String[]>();
	
	HashMap<String,StringBuilder> map = new HashMap<String,StringBuilder>();
	String linker1 = "TAGCCATCGCATTGC";
	String linker2 = "TACCTCTGAGCTGAA";
	
	public ProcessBiorad(String[] args) throws FileNotFoundException, IOException {
		if(args.length < 2)
		{
			System.err.println("Will split a biorad fastq file set (r1 and r2) into single-end fastq files for each cell, with each read annotated by the UMI and remove any non-matchd reads if arg3 is specified\n"
					+ "\nUsage: ProcessBiorad Read1.fastq.gz Read2.fastq.gz [sort?]");
			System.exit(1);
		}
		for(String barcode : barcodes)
		{	
			for(int i = 0; i < 6; i++)
			{
				char[] seq = barcode.toCharArray();
				seq[i]='A';
				closestCode.put(new String(seq), barcode);
				seq[i]='C';
				closestCode.put(new String(seq), barcode);
				seq[i]='G';
				closestCode.put(new String(seq), barcode);
				seq[i]='T';
				closestCode.put(new String(seq), barcode);
				seq[i]='N';
				closestCode.put(new String(seq), barcode);
			}
		}
		for(int i = 0; i < linker1.length();i++)
		{
			char[] seq = linker1.toCharArray();
			seq[i]='A';
			closestLinker1.add(new String(seq));
			seq[i]='C';
			closestLinker1.add(new String(seq));
			seq[i]='G';
			closestLinker1.add(new String(seq));
			seq[i]='T';
			closestLinker1.add(new String(seq));
			seq[i]='N';
			closestLinker1.add(new String(seq));
		}
		for(int i = 0; i < linker2.length();i++)
		{
			char[] seq = linker2.toCharArray();
			seq[i]='A';
			closestLinker2.add(new String(seq));
			seq[i]='C';
			closestLinker2.add(new String(seq));
			seq[i]='G';
			closestLinker2.add(new String(seq));
			seq[i]='T';
			closestLinker2.add(new String(seq));
			seq[i]='N';
			closestLinker2.add(new String(seq));
		}
		for(int i = 0; i < 3;i++)
		{
			char[] seq = "ACG".toCharArray();
			seq[i]='A';
			ed1ToACG.add(new String(seq));
			seq[i]='C';
			ed1ToACG.add(new String(seq));
			seq[i]='G';
			ed1ToACG.add(new String(seq));
			seq[i]='T';
			ed1ToACG.add(new String(seq));
			seq[i]='N';
			ed1ToACG.add(new String(seq));
		}
		for(int i = 0; i < 6;i++)
		{
			char[] seq = "GACTTT".toCharArray();
			seq[i]='A';
			ed1ToGACTTT.add(new String(seq));
			seq[i]='C';
			ed1ToGACTTT.add(new String(seq));
			seq[i]='G';
			ed1ToGACTTT.add(new String(seq));
			seq[i]='T';
			ed1ToGACTTT.add(new String(seq));
			seq[i]='N';
			ed1ToGACTTT.add(new String(seq));
		}

		Scanner ss1 = new Scanner(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]))));
		Scanner ss2 = new Scanner(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1]))));
		if(args.length > 2)
		{
			System.out.println("Reading in all reads into memory and outputing to sorted fastq.gz files which will be used for downstream processing...");
			PrintWriter rpw1 = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[0].substring(0, args[0].indexOf(".fastq.gz"))+"_sorted.fastq.gz")));
			PrintWriter rpw2 = new PrintWriter(new GZIPOutputStream(new FileOutputStream(args[1].substring(0, args[1].indexOf(".fastq.gz"))+"_sorted.fastq.gz")));
			int count = 0;
			while(ss1.hasNextLine() || ss2.hasNextLine())
			{
				String name1 = null;
				String name2 = null;
				String[] read1 = new String[4];
				String[] read2 = new String[4];
				for(int i = 0; i < 4; i++)
				{
					if(ss1.hasNextLine())
					{
						String line = ss1.nextLine();
						if(i == 0) name1 = line.split(":")[0]+line.split(":")[1]+line.split(":")[2]+line.split(":")[3]+line.split(":")[4]+line.split(":")[5]+line.split(":")[6];
						read1[i]=line;
					}
					if(ss2.hasNextLine())
					{
						String line = ss2.nextLine();
						if(i == 0) name2 = line.split(":")[0]+line.split(":")[1]+line.split(":")[2]+line.split(":")[3]+line.split(":")[4]+line.split(":")[5]+line.split(":")[6];
						read2[i]=line;
					}
				}
				if(name1 != null)
				{
					//add it
					if(reads.get(name1) == null)
						reads.put(name1, read1);
					else
					{						
						//output both at the same time!
						for(int i = 0; i < 4; i++)
						{
							rpw1.println(read1[i]);
							rpw2.println(read2[i]);
						}
						reads.remove(name1);
					}
				}
				if(name2 != null)
				{
					//add it
					if(reads.get(name2) == null)
						reads.put(name2, read2);
					else
					{						
						//output both at the same time!
						for(int i = 0; i < 4; i++)
						{
							rpw1.println(read1[i]);
							rpw2.println(read2[i]);
						}
						reads.remove(name2);
					}
				}
				count++;
				if(count % 1000000 ==0)
					System.out.print(count+" Unresolved reads: "+reads.size()+"\t\t\r");
			}
			System.out.println();
			rpw1.close();
			rpw2.close();
			ss1.close();
			ss2.close();
			reads = null;
			System.gc();
			args[0] = args[0].substring(0, args[0].indexOf(".fastq.gz"))+"_sorted.fastq.gz";			
			args[1] = args[1].substring(0, args[1].indexOf(".fastq.gz"))+"_sorted.fastq.gz";
			ss1 = new Scanner(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]))));
			ss2 = new Scanner(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1]))));
		}
		int total=0;
		int filtered=0;
		
		File dir = new File(args[1].substring(0, args[1].indexOf("_R")));
		if(dir.exists())
		{
			System.err.println("The directory "+dir.getAbsolutePath()+" already exists... will not append to existing.");
			System.exit(1);
		}
		dir.mkdir();

		while(ss1.hasNextLine())
		{
			if(!processRead(ss1,ss2))
				filtered++;
			total++;
			if(total%10000000==0)
			{
				System.out.println("Outputing reads to files... Thus far: "+(total/1000000)+" m reads");
				int count = 0;
				long avgWrite = 0;
				for(String cell: map.keySet())
				{
					PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(dir.getName()+File.separator+args[1].substring(0, args[1].indexOf("_R"))+"_Cell_"+cell+"_R1.fastq.gz",true)));
					pw.print(map.get(cell).toString());
					pw.flush();
					pw.close();
					avgWrite+=map.get(cell).length();
					map.put(cell, new StringBuilder(1024*16));
					if(count%100 ==0)
						System.out.print("Total cells so far: "+count+"\t\t\r");
					count++;
				}
				avgWrite/=count;
				System.gc();
				System.out.println("Total memory available: "+(Runtime.getRuntime().freeMemory()/(1024*1024*1024))+" Gb. Avg characters written per cell: "+avgWrite);
			}
		}
		System.out.println();
		System.out.println(String.format("Kept %d/%d or %3.3f%% Total number of cells: %d",(total-filtered),total,100.0*(total-filtered)/total, map.keySet().size()));
		int count = 0;
		for(String cell: map.keySet())
		{
			PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(dir.getName()+File.separator+args[1].substring(0, args[1].indexOf("_R"))+"_Cell_"+cell+"_R1.fastq.gz",true)));
			pw.print(map.get(cell).toString());
			pw.flush();
			pw.close();
			if(count%100 ==0)
				System.out.print("Total cells: "+count+"\t\t\r");
			count++;
		}
		System.out.println("Filtering cells with small read counts...");
		int kept = 0;
		total = 0;
		for(File f: dir.listFiles())
		{
			if(f.length() < 1000000)
			{
				f.delete();
			}
			else kept++;
			total++;
			if(total % 100 == 0)
				System.out.print("So far kept "+kept+" out of "+total+"\t\t\r");
		}
		System.out.println();
		System.out.println("Kept "+kept+" out of "+total+" cells!");
		ss1.close();
		ss2.close();
	}
	
	String[] barcodes = {"AAAGAA","AGTCTG","CCGTAA","GACTCG","GGTAGG","TCGCCT",
	"AACAGC","ATACTT","CCTCTA","GAGCTT","GGTGCT","TCGGGA",
	"AACGTG","ATAGCG","CGAAAG","GAGGCC","GTACAG","TCTAGC",
	"AAGCCA","ATATAC","CGAGCA","GAGTGA","GTCCTA","TGAATT",
	"AAGTAT","ATCCGG","CGCATA","GATCAA","GTCGGC","TGAGAC",
	"AATTGG","ATGAAG","CGGCGT","GCCAGA","GTGGTG","TGCGGT",
	"ACAAGG","ATTAGT","CGGTCC","GCCGTT","GTTAAC","TGCTAA",
	"ACCCAA","CAACCG","CGTTAT","GCGAAT","GTTTCA","TGGCAG",
	"ACCTTC","CAAGTC","CTAGGT","GCGCGG","TAAGCT","TGTGTA",
	"ACGGAC","CACCAC","CTATTA","GCTCCC","TAATAG","TGTTCG",
	"ACTGCA","CACTGT","CTCAAT","GCTGAG","TACCGA","TTAAGA",
	"AGACCC","CAGACT","CTGTGG","GCTTGT","TAGAGG","TTCGCA",
	"AGATGT","CAGGAG","CTTACG","GGACGA","TATTTC","TTCTTG",
	"AGCACG","CATAGA","CTTGAA","GGATTG","TCAGTG","TTGCTC",
	"AGGTTA","CCACGC","GAAATA","GGCCAT","TCATCA","TTGGAT",
	"AGTAAA","CCGATG","GAAGGG","GGGATC","TCCAAG","TTTGGG"};
	
	HashMap<String,String> closestCode = new HashMap<String,String>();
	HashSet<String> closestLinker1 = new HashSet<String>();
	HashSet<String> closestLinker2 = new HashSet<String>();
	HashSet<String> ed1ToACG = new HashSet<String>();
	HashSet<String> ed1ToGACTTT = new HashSet<String>();
	
	private boolean processRead(Scanner ss1, Scanner ss2) {
		String[] read1 = new String[4];
		String[] read2 = new String[4];
		for(int i = 0; i < read1.length; i++)
		{
			try{
				read1[i]=ss1.nextLine();
				read2[i]=ss2.nextLine();				
			}catch(Exception e){return false;}
		}
		//lets check to make sure that the read names match
		if(read1[0].split(":")[0].compareTo(read2[0].split(":")[0])!= 0||
				read1[0].split(":")[1].compareTo(read2[0].split(":")[1])!= 0||
				read1[0].split(":")[2].compareTo(read2[0].split(":")[2])!= 0||
				read1[0].split(":")[3].compareTo(read2[0].split(":")[3])!= 0||
				read1[0].split(":")[4].compareTo(read2[0].split(":")[4])!= 0||
				read1[0].split(":")[5].compareTo(read2[0].split(":")[5])!= 0||
				read1[0].split(":")[6].compareTo(read2[0].split(":")[6])!= 0)
		{
			System.err.println("Read1: "+read1[0]+" Read2: "+read2[0]+" don't match so skipping...");
			return false;
		}
		//lets pull the sequence from read1
		String seq = read1[1];
		//lets look for the Linker1 and Linker2 sequences
		//simplest case, 0 edit distance
		int iLink1 = seq.indexOf(linker1);
		if(iLink1 < 0)
		{
			//did not find
			for(String l: closestLinker1)
			{
				int ind = seq.indexOf(l);
				if(ind > -1)
				{
					iLink1 = ind;
					break;
				}
			}
		}
		if(iLink1 == -1) return false;  //bad read
		int iLink2 = seq.indexOf(linker2);
		if(iLink2 < 0)
		{
			//did not find
			for(String l: closestLinker2)
			{
				int ind = seq.indexOf(l);
				if(ind > -1)
				{
					iLink2 = ind;
					break;
				}
			}
		}
		if(iLink2 == -1) return false;  //bad read
		//lets filter to ensure that BC1 and BC2 are valid sizes
		if(iLink1 >= 6 && iLink2-iLink1 == 21)
		{
			String GACTTT="";
			String BC1,BC2,BC3,ACG,UMI;
			BC1=BC2=BC3=ACG=UMI="";
			try
			{
				BC1 = seq.substring(iLink1-6,iLink1);
				BC2 = seq.substring(iLink2-6,iLink2);
				BC3 = seq.substring(iLink2+15,iLink2+21);
				ACG = seq.substring(iLink2+21,iLink2+24);
				UMI = seq.substring(iLink2+24,iLink2+32);
				GACTTT= seq.substring(iLink2+32,iLink2+38);
			}catch(Exception e)
			{
				return false;
			}
			if(ed1ToACG.contains(ACG) && ed1ToGACTTT.contains(GACTTT) && UMI.indexOf('N')==-1)
			{
				if(closestCode.get(BC1) != null)
					BC1 = closestCode.get(BC1);
				else return false; //lets not try to save a weird barcode
				if(closestCode.get(BC2) != null)
					BC2 = closestCode.get(BC2);
				else return false; //lets not try to save a weird barcode
				if(closestCode.get(BC3) != null)
					BC3 = closestCode.get(BC3);
				else return false; //lets not try to save a weird barcode
				String cellBarcode = BC1+BC2+BC3;

				if(map.get(cellBarcode)==null)
				{
					map.put(cellBarcode,new StringBuilder(1024*16));
				}
				map.get(cellBarcode).append(UMI+":"+read2[0]+System.lineSeparator());
				map.get(cellBarcode).append(read2[1]+System.lineSeparator());
				map.get(cellBarcode).append(read2[2]+System.lineSeparator());
				map.get(cellBarcode).append(read2[3]+System.lineSeparator());
				return true;
			}
			else return false;
		}
		else return false;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new ProcessBiorad(args);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
