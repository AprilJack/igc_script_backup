import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;

import javax.swing.ListCellRenderer;


public class quantifyReadthrough {
	private class CDS
	{
		int start,stop;
		String id = "";
		LinkedList<Double> values = new LinkedList<Double>();
		
		public CDS(int start,int stop, String id)
		{
			this.start=start;this.stop=stop;this.id=id;
		}
		
		public void add(double value)
		{
			values.add(value);
		}
	}
	
	public quantifyReadthrough(String[] args) throws FileNotFoundException
	{
		if(args.length < 3)
		{
			System.err.println("Usage: quantifyReadthrough File.GTF window file1.bedGraph [file2.bedGraph] ...");
			System.err.println("\nWill integrate bedGraphs on each strand within the tss-window->tss, tss->tss_window, tss_window->tss_window*2"
					+ "\nfor each bedGraph and print a matrix to stdout of normalized values after filtering for overlapping genes");
			System.exit(1);
		}
		FastScanner fs = new FastScanner(args[0]);
		//load the gtf
		HashMap<String,LinkedList<CDS>> genes = new HashMap<String,LinkedList<CDS>>();
		int window=Integer.parseInt(args[1]);
		System.err.println("Loading GTF: "+args[0]+" with window= "+window);
		CDS cds = null;
		while(fs.hasMore())
		{
			String[] line = fs.getLine().split("\t");
			if(line.length < 3) continue;
			if(line[2].compareTo("exon")==0)
			{
				String chr = line[0];
			
				int start = Integer.parseInt(line[3]);
				int end = Integer.parseInt(line[4]);
				String id = line[8].substring(line[8].indexOf("transcript_id ")+14,line[8].length()-2).replaceAll("\"", "");
				if((cds==null || cds.id.compareTo(id) != 0))
				{
					if(genes.get(chr+line[6]) == null) genes.put(chr+line[6],new LinkedList<CDS>());
					cds = new CDS(start,end,id);
					genes.get(chr+line[6]).add(cds);
				}
				if(cds.id.compareTo(id) == 0)
				{
					cds.stop= Math.max(cds.stop,end);
				}
			}

		} //finish reading GTF
		for(String str: genes.keySet())
		{
			//lets filter these so they don't overlap
			boolean plus = str.charAt(str.length()-1)=='+';
			LinkedList<CDS> temp = new LinkedList<CDS>();
			CDS last = null;
			LinkedList<CDS> list = genes.get(str);
			if(!plus)
				Collections.reverse(list);
			for(CDS c: list)
			{
				//System.err.println(c.start+" "+c.stop);
				if(plus)
				{
					if(last != null && last.stop+window*2 < c.start)
					{
						//this one is ok too
						
						//temp.add(c);
						temp.add(last);
						
					}
				}
				else
				{
					if(last != null && last.start-window*2 > c.stop)
					{

						temp.add(last);
						
					}
				}
				last=c;
			}
			Collections.reverse(temp);
			genes.put(str, temp);
		}
		fs.close();
		//now lets read the bedGraph and quantify them
		for(int i = 2; i < args.length; i++)
		{
			System.err.println("Processing bedGraph: "+args[i]);
			Scanner s = new Scanner(new File(args[i]));
			LinkedList<CDS> gs = null;
			LinkedList<CDS> copy = null;
			CDS current = null;
			double sumCDS= 0;
			double sumTss = 0;
			double sumIntron =0;
			String sign = "+";
			while(s.hasNextLine())
			{
				String nextLine = s.nextLine();
			//	System.err.println(nextLine);
				if(nextLine.contains("+ strand"))
				{
					sign="+";
				}
				else if(nextLine.contains("- strand"))
				{
					sign="-";
				}
				else  //in the main body of the track
				{
					String[] line = nextLine.split("\t");
					if(line.length < 3) continue;
					String chr = line[0];
					if(genes.get(chr+sign) == null)
					{
						continue; //this is a bad chr
					}
					if(gs == null || !gs.equals(genes.get(chr+sign)))
					{
						//lets empty out the copy
						if(copy != null )
						{
							while(copy.size() > 0)
							{
								CDS remaining = copy.removeFirst();
								remaining.add(0.0);
								remaining.add(0.0);
								remaining.add(0.0);
							}
						}
						gs=genes.get(chr+sign);
						copy = new LinkedList<CDS>();
						copy.addAll(gs);
						if(current != null)
						{
							current.add(sumCDS/(Math.min(window,current.stop-current.start)));
							current.add(sumTss/window);
							current.add(sumIntron/window);
						}
						if(copy.size() > 0)
						{
							current = copy.removeFirst();
							sumCDS= 0;
							sumTss = 0;
							sumIntron =0;
						}
					}
					//now lets check if we are hitting the current cds
					int left = Integer.parseInt(line[1]);
					int right = Integer.parseInt(line[2]);
					double value=  Double.parseDouble(line[3]);
					if(sign.compareTo("+")==0)
					{
						while(left > current.stop+window*2)
						{
							//we are past the end of the windows so lets move on to the next CDS
							//first lets save what we have for this guy
							if(copy.size() > 0)
							{
								current.add(sumCDS/(Math.min(current.stop-current.start,window)));
								current.add(sumTss/window);
								current.add(sumIntron/window);
								current = copy.removeFirst();
								sumCDS=sumTss=sumIntron=0;
							}
							else break;
						}
						if(left < current.stop && right > Math.max(current.start,current.stop-window))
						{
							int length = Math.min(right,current.stop)-Math.max(Math.max(left,current.start),current.stop-window);
							sumCDS+=length*value;
						}
						else if(left < current.stop+window && right > current.stop)
						{
							int length = Math.min(right,current.stop+window)-Math.max(left,current.stop);
							sumTss+=length*value;
						}
						else if(left < current.stop+window*2 && right > current.stop+window)
						{
							int length = Math.min(right,current.stop+window*2)-Math.max(left,current.stop+window);
							sumIntron+=length*value;
						}
					}
					else  //negative strand
					{
						while(left > current.start+window)
						{
							//we are past the end of the windows so lets move on to the next CDS
							//first lets save what we have for this guy
							if(copy.size() > 0)
							{
								current.add(sumCDS/(Math.min(current.stop-current.start,window)));
								current.add(sumTss/window);
								current.add(sumIntron/window);
								current = copy.removeFirst();
								sumCDS=sumTss=sumIntron=0;
							}
							else break;
						}
						if(left < Math.min(current.start+window,current.stop) && right > current.start)
						{
							int length = Math.min(Math.min(right,current.stop),current.start+window)-Math.max(left,current.start);
							sumCDS+=length*value;
						}
						else if(left < current.start && right > current.start-window)
						{
							int length = Math.min(right,current.start)-Math.max(left,current.start-window);
							sumTss+=length*value;
						}
						else if(left < current.start-window && right > current.start-window*2)
						{
							int length = Math.min(right,current.start-window)-Math.max(left,current.start-2*window);
							sumIntron+=length*value;
						}
					}

				}
			}
			s.close();
		}
		//ok so now lets output all of the CDSs and their values for each bedGraph in one bigAss matrix
		String header ="Transcript window="+window+"\tchr\tstrand\tstart\tstop";
		for(int i =3; i < args.length; i++)
			header+="\t"+args[i]+"_CDS\t"+args[i]+"_TSS->\t"+args[i]+"_Intron->";
		System.out.println(header);
		for(String chrStr: genes.keySet())
		{
			char plus = chrStr.charAt(chrStr.length()-1);
			String chr=chrStr.substring(0,chrStr.length()-1);
			LinkedList<CDS> list = genes.get(chrStr);
			for(CDS c: list)
			{
				if(c.values.size() ==0) break;
				String output= c.id+"\t"+chr+"\t"+plus+"\t"+c.start+"\t"+c.stop;
				for(Double value: c.values)
					output+=String.format("\t%3.4f",value*100.0);
				System.out.println(output);
			}
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new quantifyReadthrough(args);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	


}
