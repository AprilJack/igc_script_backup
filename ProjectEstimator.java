import java.io.File;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;


public class ProjectEstimator {

	public ProjectEstimator(String[] args) {
		if(args.length < 1)
		{
			System.err.println("Usage: ProjectEstimator directory [directory2]...");
			System.err.println("\nWill summarize disk usage and time spent for each directory");
			System.exit(1);
		}
		System.out.println("Directory\tSize (Gb)\tDuration(days)");
		System.err.println("Directory\tSize (Gb)\tDuration(days)");
		for(String arg: args)
		{
			File ff = new File(arg);
			if(ff.isDirectory())
			{
				HashMap<String,long[]> map = new HashMap<String,long[]>();
				populateMap(new File(arg),map);
				//now lets characterize the files
				LinkedList<Long> dates = new LinkedList<Long>();
				
				double totalSize = 0; //in gigs
				for(String f: map.keySet())
				{
					totalSize+=map.get(f)[0]/(1024.0*1024.0*1024.0);
					dates.add(map.get(f)[1]);
				}
				Collections.sort(dates);
				if(dates.size() > 1)
				{
					long last = dates.getFirst();
					double duration = 0;
					for(long date: dates)
					{
						double diff = (date-last)/(24*3600000.0);
						if(diff< 10)
							duration+=diff;
						last = date;
					}
					if(duration > 0.01 && totalSize > 1)
					{
						System.out.println(String.format("%s\t%3.3f\t%3.3f",arg,totalSize,duration));
						System.err.println(String.format("%s\t%3.3f\t%3.3f",arg,totalSize,duration));
					}

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
		new ProjectEstimator(args);
	}

}
