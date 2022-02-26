import java.util.ArrayList;
import java.util.HashMap;


public class bedCoverageBias {

	public bedCoverageBias(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: bedCoverageBias bedFileWithLocationsOfInterest(genes) samFile");
			System.err.println("Outputs the average % coverage difference from %50 of the regions in the bed file");
			System.err.println("Only samples a million reads from the sam file to keep it fast");
			System.exit(1);
		}
		//read in the entire bed file
		HashMap<String,ArrayList<long[]>> chrs = new HashMap<String,ArrayList<long[]>>();
		FastScanner fs = new FastScanner(args[0]);
		HashMap<String,double[]> coverage = new HashMap<String,double[]>();
		while(fs.hasMore())
		{
			String line = fs.getLine();
			if(!line.startsWith("#"))
			{
				String[] split = line.split("\t");
				if(split.length > 2)
				{
					if(chrs.get(split[0]) == null) chrs.put(split[0], new ArrayList<long[]>());
					long[] feature = new long[3];
					feature[0] = Long.parseLong(split[1]);
					feature[1] = Long.parseLong(split[2]);
					if(split[5].compareTo("+")==0)
						feature[2] = 1;
					else
						feature[2] = 0;
					int index = 0;
					while(index < chrs.get(split[0]).size())
					{
						if(chrs.get(split[0]).get(index)[0] >feature[0])
							break;
						index++;
					}
					chrs.get(split[0]).add(index,feature);
				}
			}
		}
		int count = 0;
		//reads in the sam file
		FastScanner s = new FastScanner(args[1]);
		while(s.hasMore())
		{
			String line = s.getLine();
			if(line.startsWith("@"))
			{
				//header we can ignore...
			}
			else
			{
				if(count == 1000000)
					break;
				count++;
				String[] split = line.split("\t");
				int flag = Integer.parseInt(split[1]);
				boolean forward = true;
				if(Integer.toBinaryString(flag).charAt(4)==1)
					forward = false;
				String chr = split[2];
				int start = Integer.parseInt(split[3]);
				String cigar = split[5];
				String[] puffs = cigar.split("[MNSX=HPD]");
				int stop = start;

				for(int i = 0; i < puffs.length; i++)
				{
					try{
						stop+=Integer.parseInt(puffs[i]);
					}catch(Exception e){}
				}
				double average = start+(stop-start)/2.0;
				ArrayList<long[]> features = chrs.get(chr);
				if(features != null)
				{
					for(int i = 0; i < features.size(); i++)
					{
						if(average >= features.get(i)[0] && average <= features.get(i)[1])
						{
							if((forward && features.get(i)[2] ==1)||(!forward && features.get(i)[2]==0))
							{
								//same strand
								double bias = ((stop+start)/2.0-features.get(i)[0])/(features.get(i)[1]-features.get(i)[0]);
								if(!forward)
									bias = 1-bias;
								String region = features.get(i)[0]+","+features.get(i)[1];
								if(coverage.get(region) == null)
									coverage.put(region, new double[2]);
								coverage.get(region)[0]+=bias;
								coverage.get(region)[1]+=1;
							}
						}
					}
				}
			}
		}
		//now lets output the average and standard deviation of the biases across these features
		double avg = 0;
		count = 0;
		for(String str: coverage.keySet())
		{
			if(coverage.get(str)[1] >= 10)
			{
				avg+=coverage.get(str)[0]/coverage.get(str)[1];
				count++;
			}
		}
		avg/=count;
		double sd = 0;
		for(String str: coverage.keySet())
		{
			if(coverage.get(str)[1] >= 10)
			{
				double val =coverage.get(str)[0]/coverage.get(str)[1];
				sd+= (val-avg)*(val-avg);
			}
		}
		sd = Math.sqrt(sd/count);
		System.out.println(args[0]+"\t"+args[1]+"\t"+avg+"\t"+sd+"\t"+(avg-sd/Math.sqrt(count))+"\t"+(avg+sd/Math.sqrt(count))+"\t"+count);
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new bedCoverageBias(args);
	}

}
