import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;


public class plotDomains {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 1)
		{
			System.err.println("Usage: plotDomains HiCMap.txt [skipChr1] [skipSChr2] [...] > HiCMap.bedGraph");
			System.exit(1);
		}
		FastScanner s1;
		try {
			System.err.println("Loading "+args[0]);
			s1 = new FastScanner(args[0]);
			String color = ((int)(Math.random()*200))+","+((int)(Math.random()*200))+","+((int)(Math.random()*200));
			HashSet<String> skipChr = new HashSet<String>();
			if(args.length > 1)
			{
				for(int i = 1; i < args.length; i++)
					skipChr.add(args[i]);
			}
			int j = 0;
			s1.getLine();
			HashMap<String,double[]> histograms = new HashMap<String,double[]>();
			HashMap<String,Integer> chrStarts = new HashMap<String,Integer>();
			//HashMap<String,double[]> rows = new HashMap<String,double[]>();
			HashMap<String,ArrayList<Integer>> positions = new HashMap<String,ArrayList<Integer>>();
			double distances = 0.1; 
			double totalCounts = 0.1;
			String currentChr = null;
			int chrLength =0;
			while (s1.hasMore())
			{
				String[] line1 = s1.getLine().split("\t");
				String chr = line1[0].split("-")[0];
				if(currentChr == null)
				{
					currentChr = chr;
					chrStarts.put(chr,0);
				}
				else if(currentChr.compareTo(chr)!=0)
				{
					//prep chr
					histograms.put(currentChr, new double[chrLength]);
					currentChr = chr;
					chrLength=0;
					chrStarts.put(chr,j);
				}
				chrLength++;
				j++;
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}
			System.err.println();
			histograms.put(currentChr, new double[chrLength]);
			s1.close();
			j=0;
			String lastChr = null;
			s1=new FastScanner(args[0]);
			s1.getLine();
			while (s1.hasMore())
			{
				String[] line1 = s1.getLine().split("\t");
				if(line1.length < 3) continue;
				String chr = line1[0].split("-")[0];
				if(chr.compareTo("chr1")==0)
				{
					System.err.print("");
				}
				if(lastChr == null) lastChr=chr;
				if(lastChr.compareTo(chr) != 0)
				{
					totalCounts/=(0.1+Math.pow(line1.length-2,2));
					for(int i = 0; i < histograms.get(chr).length; i++)
					{
						histograms.get(lastChr)[i]/=((totalCounts)*distances);
					}
					totalCounts=0.1;
					distances=0.1;
				}
				int position = Integer.parseInt(line1[0].split("-")[1]); //but this is genomic position
				int chrPos = 0;
				for(int i = chrStarts.get(chr)+2; i < chrStarts.get(chr)+2+histograms.get(chr).length; i++)
				{
					double value = Double.parseDouble(line1[i]);
					totalCounts+=value;
					//if(value > 1) value = 1;
					if(Double.isInfinite(value))
						value = 0; //not sure what to do here (the input is messed up)
					histograms.get(chr)[i-2-chrStarts.get(chr)]+=value*(1+Math.abs((i-2)-j));
					distances+=((double)(1+Math.abs((i-2)-j)))/(line1.length-2);
				}
				if(positions.get(chr) == null)
				{
					positions.put(chr, new ArrayList<Integer>());
					chrPos = position; 
				}
				positions.get(chr).add(position-chrPos);

				j++;			
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}

			//now lets rescale the histogram
			System.err.println();
			int window=100;
			System.out.println("track type=bedGraph name=\"Domains."+args[0]+"\" yLineMark=0.0 alwaysZero=on visibility=full autoScale=on color="+color);

			for(String chr : histograms.keySet())
			{
				if(skipChr.contains(chr))
				{
					System.err.println("Skipping "+chr);
					continue;
				}
				else
					System.err.println("Processing "+chr);
				double[] output = new double[histograms.get(chr).length];
				double last = 0;
				int prevPos=1;
				for(int i = 0; i < histograms.get(chr).length; i++)
				{
					int count = 0;
					double val = 0;
					
					for(int w = -window; w <= 0; w++)
					{
						try{
							val+=histograms.get(chr)[i+w];
							count++;
						}catch(Exception e){count++;}
					}
					
					double average=val/count;
					output[i]=average-last;
					if(i > 0 && i < positions.get(chr).size())
					{
						if(Double.isNaN(average-last))
						{
							average = 0;
							last = 0;
						}
						System.out.println(String.format("%s\t%d\t%d\t%3.6f",chr,prevPos,positions.get(chr).get(i),average-last));
						prevPos=positions.get(chr).get(i);
					}
					last = average;
						
				}
				
			//		System.err.println(rows.get(chr)[i]);
			}
			
			//lets average everything out
			/*
			for(String chr : rows.keySet())
			{
//				for(double threshold = 0.1)
				ArrayList<ArrayList<Double>> values = rows.get(chr);
				ArrayList<ArrayList<Double>> clusterAverages = new ArrayList<ArrayList<Double>>();
				System.err.println("Finished loading data");
				for(int c = 0; c < 2; c++)
				{
					ArrayList<Double> temp = new ArrayList<Double>();
					for(int i = 0; i < values.size(); i++)
					{
						temp.add(0.0);
					}
					clusterAverages.add(temp);
				}
				int id = 0;
				//double[] pearsons = new double[values.size()];
				int lastSwitch = 0;
				for(int i = 0; i < values.size()-1; i++)
				{
					ArrayList<Double> current = values.get(i);
					double avgPearson = 0;
					int count = 0;
					try{
						for(int n = 1; n <= 3; n++)
						{
							ArrayList<Double> next = values.get(i+n);
							double xbar = 0;
							for(double val: current)
							{
								xbar+=val;
							}
							double ybar = 0;
							for(double val: next)
							{
								ybar+=val;
							}
							xbar/=current.size();
							ybar/=next.size();
							double top = 0;
							double sqrt1= 0;
							double sqrt2 =0;
							for(int x = 0; x < current.size(); x++)
							{
								top+=(current.get(x)-xbar)*(next.get(x)-ybar);
								sqrt1 += Math.pow(current.get(x)-xbar,2);
								sqrt2 += Math.pow(next.get(x)-ybar,2);
							}
							double pearson = top/(Math.sqrt(sqrt1)*Math.sqrt(sqrt2));
							avgPearson+=pearson;
							count++;
						}
					}catch(Exception e){}
					avgPearson/=count;
					//pearsons[i+1]=pearson;
					if(avgPearson < 0.30 && i-lastSwitch > 3) // this match sucks
					{
						id = Math.abs(id-1);
						lastSwitch =i;
					}

					for(int k =0; k < clusterAverages.get(id).size(); k++)
					{
						clusterAverages.get(id).set(k,clusterAverages.get(id).get(k)+current.get(k));
					}
				
					if(id == 0)
						System.out.println(chr+"\t"+1+"\t"+positions.get(chr).get(i)+"\t-1");
					else
						System.out.println(chr+"\t"+1+"\t"+positions.get(chr).get(i)+"\t1");
				
					//if(i%(values.size()/100)==0)
					//	System.err.print(".");
				}
				
				//now lets compute the score for each row again
				System.out.println("track name=\"sdDist."+args[0]+"\" yLineMark=\"0.0\" alwaysZero=on maxHeightPixels=100:75:11 visibility=full viewLimits=-1:1 autoScane = off type=bedGraph");
				int prevPos=1;
				for(int i = 1; i < values.size(); i++)
				{
					ArrayList<Double> current = values.get(i);
					//now lets compute the pearson toward cluster1
					double[] pearsons = new double[clusterAverages.size()];
					for(int c = 0; c < clusterAverages.size();c++)
					{
						ArrayList<Double> next = clusterAverages.get(c);
						double xbar = 0;
						for(double val: current)
						{
							xbar+=val;
						}
						double ybar = 0;
						for(double val: next)
						{
							ybar+=val;
						}
						xbar/=current.size();
						ybar/=next.size();
						double top = 0;
						double sqrt1= 0;
						double sqrt2 =0;
						for(int x = 0; x < current.size(); x++)
						{
							top+=(current.get(x)-xbar)*(next.get(x)-ybar);
							sqrt1 += Math.pow(current.get(x)-xbar,2);
							sqrt2 += Math.pow(next.get(x)-ybar,2);
						}
						double pearson = top/(Math.sqrt(sqrt1)*Math.sqrt(sqrt2));
						pearsons[c]=pearson;
					}
					
					double score = pearsons[0]-pearsons[1];
					System.out.println(chr+"\t"+prevPos+"\t"+positions.get(chr).get(i)+"\t"+score);
					prevPos=positions.get(chr).get(i);
					//System.err.println(i+"\t"+clusterAverages.get(0).get(i)+"\t"+clusterAverages.get(1).get(i));
				}
			}
			*/
		}
		catch(Exception e)
		{
		// TODO Auto-generated catch block
		
			e.printStackTrace();
		}
		
	}

}
