import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;


public class ProjectEstimatorTime {

	private class DateSize {
		long date=0;
		double size=0;
		public DateSize(long date, double size)
		{
			this.date=date; this.size=size;
		}

	}
	
	public ProjectEstimatorTime(String[] args) {
		if(args.length < 1)
		{
			System.err.println("Usage: ProjectEstimator directory [directory2]...");
			System.err.println("\nWill summarize disk usage and time spent for each directory. ");
			System.exit(1);
		}
		System.out.println("Directory\tSize (Gb)\tDuration(days)");
		for(String arg: args)
		{
			File ff = new File(arg);
			if(ff.isDirectory())
			{
				System.err.print(arg);
				HashMap<String,long[]> map = new HashMap<String,long[]>();
				populateMap(new File(arg),map);
				//now lets characterize the files
				LinkedList<DateSize> dates = new LinkedList<DateSize>();
				double totalSize = 0; //in gigs
				for(String f: map.keySet())
				{
					double size=map.get(f)[0]/(1024.0*1024.0*1024.0);
					dates.add(new DateSize(map.get(f)[1],size));
					totalSize+=size;
				}
				Collections.sort(dates, new Comparator<DateSize>(){

					@Override
					public int compare(DateSize o1, DateSize o2) {
						// TODO Auto-generated method stub
						return (int) (o1.date/1000.0-o2.date/1000.0);
					}
					
				});
				if(dates.size() > 1)
				{
					DateSize last = dates.getFirst();
					
					double duration = 0;
					for(DateSize date: dates)
					{
						double diff = (date.date-last.date)/(24*3600000.0);
						//System.err.println("date:"+date.date);
						if(diff< 5)
							duration+=diff;
						last = date;
					}
					if(duration > 0.01 && totalSize > 0.01)
					{
						System.out.println(String.format("%s\t%3.3f\t%3.3f",arg,totalSize,duration));
					}
					int day = 0;
					long now = System.currentTimeMillis();
					Collections.reverse(dates);
					double cumulative = totalSize;
					for(DateSize date: dates)
					{
						cumulative-=date.size;
						int diff = (int)((now-date.date)/(24*3600000.0));
						//System.err.println("Diff:"+diff);
						while(diff >= day)
						{
							System.err.print(String.format(",%3.2f",cumulative));
							day++;
						}
					}
					System.err.println(",0");
				}
				map.clear();
			}
		}
		
	}

	private void populateMap(File f, HashMap<String,long[]> map)
	{
		try{
			if(f.isDirectory())
			{
				for(File ff: f.listFiles())
				{
					if(map.get(ff.getAbsolutePath())==null)
						populateMap(ff,map);
				}
			}
			else
			{
	
				long[] vals = {f.length(),f.lastModified()};
				if(vals[1]!=0)
					map.put(f.getAbsolutePath(),vals);
			}
		}catch(Exception e){}
	}


	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ProjectEstimatorTime(args);
	}

}
