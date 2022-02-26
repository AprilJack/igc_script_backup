import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;


public class addColumns {

	HashMap<String,ArrayList<String>> map = new HashMap<String,ArrayList<String>>();
	LinkedList<String> order = new LinkedList<String>();
	public addColumns(String[] args) throws FileNotFoundException, IOException {
		//lets combine all of the files in the args one by one into one file and print it out
		if(args.length == 0)
		{
			System.err.println("Usage: addColumns file1.tab file2.tab.gz [file3.tab.gz] ... > matrix_sorted_by_first_column");
			System.err.println("Usage2: addColumns directory_with_tab_delimited_files_only > matrix_sorted_by_first_column");
		}
		if(args.length == 1)
		{
			//we are loading in a dir
			File dir =new File(args[0]);
			loadColumns(dir.list(), args[0]);
		}
		else
			loadColumns(args,"");
		//now lets output the matrix as one block
		for(String header: order)
		{
			String out = header;
			for(String value: map.get(out))
				out+="\t"+value;
			System.out.println(out);
		}
	}
	
	private void loadColumns(String[] args, String dir) throws FileNotFoundException, IOException
	{
		for(String input: args)
		{
			if(dir.length() != 0)
				input = dir+File.separator+input;
			Scanner s = null;
			if(input.endsWith(".gz"))
				s = new Scanner(new GZIPInputStream(new FileInputStream(input)));
			else
				s = new Scanner(new File(input));
			//now lets read it in
			while(s.hasNextLine())
			{
				String[] split = s.nextLine().split("\t");
				if(map.get(split[0])==null){
					map.put(split[0],new ArrayList<String>());
					order.add(split[0]);
				}
				for(int i = 1; i < split.length; i++)
					map.get(split[0]).add(split[i]);
			}
			s.close();
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			new addColumns(args);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
