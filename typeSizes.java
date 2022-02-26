import java.io.File;
import java.util.HashMap;



public class typeSizes {

	HashMap<String,Integer> extensions = new HashMap<String,Integer>();
	HashMap<String,Integer> textensions = new HashMap<String,Integer>();
	public typeSizes(String[] dirs)
	{
		HashMap<String,Integer> totals = new HashMap<String,Integer>();
		for(String d: dirs)
		{
			if((new File(d)).exists())
			{
				System.out.println(d);
				recursiveSearch(new File(d));
				//lets output the top extensions by size
				int total = 0;
				for(int i = 0; i < 10; i++)
				{
					int best = 0;
					String bestExt = "";
					for(String ext: extensions.keySet())
					{
						if(extensions.get(ext) > best)
						{
							best = extensions.get(ext);
							bestExt = ext;
						}
					}
					if(best > 1024*1024)
						System.out.println(String.format("Extension: %s total size: %3.3f G",bestExt,best/(1024.0*1024.0)));
					else if(best > 1024)
						System.out.println(String.format("Extension: %s total size: %3.3f M",bestExt,best/(1024.0)));
					else
						System.out.println("Extension: "+bestExt+" total file size: "+best+" K");
					extensions.remove(bestExt);
					total+=best;
				}
				while(extensions.size() > 0)
				{
					total+= extensions.remove(extensions.keySet().iterator().next());
				}
				totals.put(d,total);
				total = 0;

			}else{System.err.println("Error: "+d+" doesn't seem to exist");}
		}
		System.out.println("Top extensions:");
		for(int i = 0; i < 10; i++)
		{
			int best = 0;
			String bestExt = "";
			for(String ext: textensions.keySet())
			{
				if(textensions.get(ext) > best)
				{
					best = textensions.get(ext);
					bestExt = ext;
				}
			}
			System.out.println(String.format("Extension:\t%s\ttotal size:\t%3.3f\tG",bestExt,best/(1024.0*1024.0)));
			textensions.remove(bestExt);
		}
		System.out.println("Totals:");
		for(String d: totals.keySet())
		{
			int total = totals.get(d);
			System.out.println(String.format("%s\t%3.3f\tG",d,total/(1024.0*1024.0)));
		}

	}
	
	private void recursiveSearch(File file) {
		if(file.isDirectory())
		{
			for(File f: file.listFiles())
				recursiveSearch(f);
		}
		else
		{
			if(file.getName().contains("."))
			{
				String extension = file.getName().substring(file.getName().indexOf('.'),file.getName().length());
				if(extensions.get(extension) == null) extensions.put(extension,(int)(file.length()/(1024)));
				else extensions.put(extension, extensions.get(extension)+(int)(file.length()/(1024)));
				if(textensions.get(extension) == null) textensions.put(extension,(int)(file.length()/(1024)));
				else textensions.put(extension, textensions.get(extension)+(int)(file.length()/(1024)));
			}
		}
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new typeSizes(args);
	}

}
