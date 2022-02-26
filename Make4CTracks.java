import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;


public class Make4CTracks {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 1)
		{
			System.err.println("Usage: Make4CTracks HiCMap.txt > HiCMap.bedGraph");
			System.exit(1);
		}
		Scanner s1;
		try {
			s1 = new Scanner(new File(args[0]));
			int j = 0;
			s1.nextLine();
			double[] histogram = null;
			HashMap<String,double[]> rows = new HashMap<String,double[]>();
			HashMap<String,ArrayList<Integer>> positions = new HashMap<String,ArrayList<Integer>>();
			while (s1.hasNextLine())
			{
				String[] line1 = s1.nextLine().split("\t");
				String chr = line1[0].split("-")[0];
				int position = Integer.parseInt(line1[0].split("-")[1]);
				if(rows.get(chr) == null)
				{
					histogram = new double[line1.length-2];
					rows.put(chr, histogram);
				}
				for(int i = 2; i < line1.length; i++)
				{
					double value = Double.parseDouble(line1[i]);

					//if(value > 1) value = 1;
					histogram[i-2]+=value*(1+Math.abs((i-2)-j));
				}
				if(positions.get(chr) == null) positions.put(chr, new ArrayList<Integer>());
				positions.get(chr).add(position);
				j++;			
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}
			System.err.println();
			System.out.println("track name=\"sdDist."+args[0]+"\" yLineMark=\"0.0\" alwaysZero=on maxHeightPixels=100:75:11 visibility=full viewLimits=-1:1 autoScane = off type=bedGraph");

			int window=100;
			int prevPos=1;
			for(String chr : rows.keySet())
			{
				double last = 0;
				for(int i = 0; i < rows.get(chr).length; i++)
				{
					int count = 0;
					double val = 0;
					
					for(int w = -window; w <= 0; w++)
					{
						try{
							val+=rows.get(chr)[i+w];
							count++;
						}catch(Exception e){}
					}
					
					double average=val/count;
					if(i > 0)
						System.out.println(String.format("%s\t%d\t%d\t%3.2f",chr,prevPos,positions.get(chr).get(i),average-last));
					last = average;
					prevPos=positions.get(chr).get(i);
						
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
