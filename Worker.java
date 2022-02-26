package polyATrimmer;

import java.util.LinkedList;

public class Worker extends Thread
{
	LinkedList<String> rs;
	LinkedList<String> trimmed = new LinkedList<String>();
	String[] start;
	String[] end;
	int misMatch = 2;
	int trim,removed,total;
	double qcutoff = 0.0;
	public Worker(LinkedList<String> linkedList, String[] start, String[] end, int misMatch, double qcutoff)
	{
		this.rs = linkedList;
		this.start = start;
		this.end = end;
		this.misMatch = misMatch;
		this.qcutoff=qcutoff;
		this.start();
	}
	public void run()
	{
		while(rs.size() > 0)
		{
			String[] split = rs.removeFirst().split("}}");
			String name = split[0];
			String seq = split[1];
			String sep = split[2];
			String qual = split[3];
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
				trimmed.add(name+"}}"+seq+"}}"+sep+"}}"+qual);
				if(gotTrimmed)
					trim++;
			}
			else{
				trim++;removed++;
				trimmed.add("TRIMMED!");
			}
			total++;
			//if(((total-1)*100)/(rs.size()+trimmed.size()) < ((total)*100)/(rs.size()+trimmed.size())) System.err.println("Thread:"+toString()+" finished "+((total)*100)/(rs.size()+trimmed.size())+"%");
		}
		System.err.println("Thread: "+toString()+" completed trimming!");
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
}

