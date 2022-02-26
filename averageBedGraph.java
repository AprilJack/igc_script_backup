
public class averageBedGraph {

	public averageBedGraph(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: averageBedGraph file.bedGraph resolution.\nOutputs new averaged bedGraph file to stdout");
			System.exit(1);
		}
		FastScanner fs = new FastScanner(args[0]);
		int resolution = Integer.parseInt(args[1]);
		int currentRegion = 0;
		int lastStop = 0;
		double regionValue = 0;
		String currentChr = null;
		while(fs.hasMore())
		{
			String line = fs.getLine();
			try{
				String[] split = line.split("\t");
				if(split.length ==4)
				{
					String chr = split[0];
					int start = Integer.parseInt(split[1]);
					int stop = Integer.parseInt(split[2]);
					double value = Double.parseDouble(split[3]);
					if(currentChr == null || currentChr.compareTo(chr) != 0)
					{
						if(currentChr != null)
						{
							//lets output the last region we were working on
							if(lastStop-currentRegion*resolution > 0)
								System.out.println(currentChr+"\t"+currentRegion*resolution+"\t"+lastStop+"\t"+regionValue/(lastStop-currentRegion*resolution));
							lastStop =0;
							currentRegion = 0;
							regionValue = 0;
						}
						currentChr = chr; 
					}
					if(stop > (currentRegion+1)*resolution)
					{
						//this bedGraph spans the boundary 
						regionValue+=((currentRegion+1)*resolution-start)*value;
						System.out.println(currentChr+"\t"+currentRegion*resolution+"\t"+(currentRegion+1)*resolution+"\t"+regionValue/resolution);
						currentRegion++;
						regionValue = (stop-currentRegion*resolution)*value;
					}
					else
					{
						regionValue+=(stop-start)*value;
					}
					lastStop = stop;
				}
				else System.out.println(line);
			}catch(Exception e){System.out.println(line);}
		}
		
	}
	
	public static void main(String[] args) {
		new averageBedGraph(args);
	}

}
