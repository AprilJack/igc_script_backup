package simpleTrimmer;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class ReadTrimmer {

	public static void main(String[] args) {
		if(args.length < 1)
		{
			printUsage();
			System.exit(1);
		}
		// TODO Auto-generated method stub
		String[] end = {"AAAAAAAAAAAAAAAAAAAA"};
		String[] start = {};
		int misMatch = 2;
		int quality = 66;  //top 33% kept...
		int minq=9999;
		int maxq= 0;
		int trimmed = 0;
		int removed = 0;
		int total = 0;
		try {
			InputStream is = null;
			if(args[0].endsWith(".gz"))
				is = new GZIPInputStream(new FileInputStream(args[0]),1024*1024);
			else
				is = new FileInputStream(args[0]);
			Scanner s = new Scanner(is);
			Scanner s2 = null;
			if(args.length  > 1)
			{
				int ptr = 1;
				while(args.length < ptr)
				{
					if(args[ptr].compareTo("-3")==0 && args[ptr+1].compareTo("-")!= 0)
					{
						end =args[ptr+1].split(",");
						ptr+=2;
					}
					else if(args[ptr].compareTo("-5")==0 && args[ptr+1].compareTo("-")!= 0)
					{
						ptr+=2;
						start =args[ptr+1].split(",");
					}
					else if(args[ptr].compareTo("-m")==0)
					{
						ptr+=1;
						misMatch =Math.max(0,Integer.parseInt(args[ptr+1]));
					}
					else if(args[ptr].compareTo("-q")==0)
					{
						ptr+=1;
						quality =Math.max(0,Integer.parseInt(args[ptr+1]));
					}
					else
					{
						s2 = new Scanner(new File(args[ptr]));
						ptr++;
					}
				}
			}
			for(String str: end)
				System.err.println("Trimming: "+str+" from trailing end");
			for(String str: start)
				System.err.println("Trimming: "+str+" from beginning");


			if(s2 == null)
			{
				int samples = 0;
				while(s.hasNextLine() && samples < 100000)
				{
					String header = s.nextLine();
					String seq = s.nextLine();
					String sep = s.nextLine();
					String qual = s.nextLine();
					for(Byte b : qual.getBytes())
					{
						if(b < minq) minq = b;
						else if(b > maxq ) maxq = b;
					}
					samples++;
				}
				if(args[0].endsWith(".gz"))
					is = new GZIPInputStream(new FileInputStream(args[0]),1024*1024);
				else
					is = new FileInputStream(args[0]);
				s = new Scanner(is);
				double qcutoff =minq+(((maxq-minq)*(quality))/100.0);
				System.err.println("Min Quality was "+((char)minq)+ " and Max Quality was "+((char)maxq));
				System.err.println("Trimming reads with quality < "+((char)((int)qcutoff))+" from trailing end.");
				GZIPOutputStream out = new GZIPOutputStream(System.out,64*1024*1024);

				StringBuilder sb = new StringBuilder(64*1024*1025);
				while(s.hasNextLine())
				{
					String header = s.nextLine();
					String seq = s.nextLine();
					String sep = s.nextLine();
					String qual = s.nextLine();
					boolean gotTrimmed = false;
					for(int i = qual.length()-1; i >= 0; i--)
					{
						if(qual.charAt(i) < qcutoff)
						{
							//still a bad nt so keep trimming
						}
						else
						{
							//quality is finally higher so lets cut off til here
							seq = seq.substring(0, i+1);
							qual = qual.substring(0,i+1);
							gotTrimmed = true;
							break;
						}
						if(i == 0)
						{
							//quality is finally higher so lets cut off til here
							seq = "";
							qual = "";
							gotTrimmed = true;
						}
					}
					if(seq.length() >= 18)
					{
						for(String str: end)
						{
							int index = matchEnd(seq,str,misMatch);
							if(index > 0)
							{
								gotTrimmed =true;
								seq = seq.substring(0, seq.length()-index);
								qual = qual.substring(0,qual.length()-index);
							}
							if(seq.length() < 18)
								break;
						}
						if(seq.length() >= 18)
						{
							for(String str: start)
							{
								int index = matchStart(seq,str,misMatch);
								if(index > 0)
								{
									gotTrimmed =true;
									seq = seq.substring(index, seq.length());
									qual = qual.substring(index, qual.length());
								}
								if(seq.length() < 18)
									break;
							}
						}
					}
					if(seq.length() >= 18)
					{
						sb.append(header);sb.append(System.lineSeparator());
						sb.append(seq);sb.append(System.lineSeparator());
						sb.append(sep);sb.append(System.lineSeparator());
						sb.append(qual);sb.append(System.lineSeparator());
						if(sb.length()> 64*1024*1024)
						{
							out.write(sb.toString().getBytes());
							sb = new StringBuilder(64*1024*1025);
							System.err.println(String.format("Trimmed %d of %d total. Removed %d short reads or %3.2f%%",trimmed,total,removed,((removed*100.0)/total)));
						}
						if(gotTrimmed)
							trimmed++;
					}
					else{
						trimmed++;removed++;
					}
					total++;
				}
				s.close();

				out.write(sb.toString().getBytes());
				out.flush();
				out.close();
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block

			printUsage();
			e.printStackTrace();
		}
		System.err.println("STATS: Trimmed "+trimmed+" of which "+removed+" reads were removed out of "+total);
	}

	private static int matchEnd(String read, String seq, int misMatch) {
		int[] mismatches = new int[Math.min(seq.length(),read.length())];
		for(int i = 0; i < seq.length() && i < read.length()-1; i++)
		{
			for (int j = 0; j <= i; j++)
				if(read.charAt(read.length()-1-i+j)!=seq.charAt(i-j))
					mismatches[i]++;
		}
		for(int i = mismatches.length-1; i > misMatch*2; i--)
			if(mismatches[i] <= misMatch)
				return i;
		return 0;
	}
	
	private static int matchStart(String read, String seq,int misMatch) {
		int[] mismatches = new int[Math.min(seq.length(),read.length())];
		for(int i = 0; i < seq.length() && i < read.length()-1; i++)
		{
			for (int j = 0; j <= i; j++)
				if(read.charAt(j)!=seq.charAt(seq.length()-1-i+j))
					mismatches[i]++;
		}
		for(int i = mismatches.length-1; i >misMatch*2; i--)
			if(mismatches[i] <= misMatch)
				return i;
		return 0;
	}

	private static void printUsage() {
		System.err.println("Usage: Trimmer fastq [-3 SEQ1,SEQ2,... (default: polyA)] [-5 SEQ1,SEQ2,... default empty] [-m mismatch default 2] [-q quality threshold (default = 66% of max)]");
		System.err.println("Currently only accepts 1 fastq or fastq.gz file and outputs gzipped filtered reads to stdout.");
		System.err.println("Trims seq by sliding seq from 3' or 5' end of read until mismatch > set default (default 2), removing nts as it goes.");
		System.err.println("Providing a \"-\" specifies that you don't want to trim from that end");
	}
	
	

}

