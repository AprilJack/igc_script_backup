import java.util.LinkedList;


public class findBoundaries {

	
	public findBoundaries(String[] args) {
		if(args.length < 1)
		{
			System.err.println("Usage:findBoundaries bedGraph [res=10000] [thresold=0] > peaks.txt");
			System.exit(1);
		}
		FastScanner fs = new FastScanner(args[0]);
		double threshold = 0;
		int res = 10000;
		if(args.length > 1)
			res = Integer.parseInt(args[1]);
		if(args.length > 2)
			threshold = Double.parseDouble(args[2]);
		String chr = "";
		int peak = 0;
		int lastPeakUp = -1000000000;
		int lastPeakDown = -100000000;
		System.out.println("#Boundaries of "+args[0]+" with res="+res+" threshold="+threshold);
		System.out.println("#PeakID\tchr\tleft\tright\t+\tChange in bedGraph");
		LinkedList<Double> values = new LinkedList<Double>();
		double avg = 0;
		while(fs.hasMore())
		{
			String line = fs.getLine();
			String[] split = line.split("\t");
			if(split.length > 3)
			{
				double value = Double.parseDouble(split[3]);
				int left = Integer.parseInt(split[1]);
				int right = Integer.parseInt(split[2]);
				if(split[0].compareTo(chr) != 0)
				{
					chr = split[0];
					//restart
					values = new LinkedList<Double>();
					avg=0;
					peak = 0;
					lastPeakUp = -10000000;
					lastPeakDown = -100000000;
				}
				for(int i = left; i <= right; i++)
				{
					double oldAvg = avg;
					values.add(value);
					avg=(avg*(values.size()-1)+value)/values.size();
					if(values.size() > res)
					{
						double first =values.removeFirst();
						avg=(avg*(values.size()+1)-first)/values.size();
					}
					if(((avg > threshold && oldAvg < threshold)) && values.size()> 100 && Math.abs(i-Math.max(lastPeakUp,lastPeakDown)) > res/2)
					{
						peak++;
						//this is a boundary
						
						System.out.println(chr+"_"+peak+"\t"+chr+"\t"+(int)Math.max(1,i-res/(2)-100)+"\t"+(int)Math.max(1,i-res/(2)+100)+"\t+\t"+res*(avg-oldAvg));
						lastPeakUp = i;

					}
					else if(((avg < threshold && oldAvg > threshold)) && values.size()> 100 && Math.abs(i-Math.max(lastPeakUp,lastPeakDown)) > res/2)
					{
						peak++;
						//this is a boundary
						//double rate = Math.max(0.1,Math.abs(res*(avg-oldAvg)));

						System.out.println(chr+"_"+peak+"\t"+chr+"\t"+(int)Math.max(1,i-res/(2)-100)+"\t"+(int)Math.max(1,i-res/(2)+100)+"\t-\t"+res*(avg-oldAvg));
						lastPeakDown = i;

					}
				}
			}
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new findBoundaries(args);
	}

}
