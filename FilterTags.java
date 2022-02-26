import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class FilterTags {

	public FilterTags(String[] args) {
		if(args.length < 2)
		{
			System.err.println("Usage: FilterTags tags.tsv toChr > filtered.tsv");
			System.err.println("Ex 1: FilterTags chr1.tsv chr2 > chr1-to-chr2-only.tsv");
			System.err.println("Ex 2: FilterTags crh1.tsv -chr2 > chr1-to-allbutchr2.tsv");
			System.exit(1);
		}
		String file = args[0];
		String chr1 = args[1];
		boolean invert = false;
		if(chr1.startsWith("-"))
		{
			invert = true;
			chr1 = chr1.substring(1);
		}
		try {
			Scanner s = new Scanner(new File(file));
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\t");
				boolean matches1 = split[6].compareTo(chr1)==0;
				if(invert) matches1 = !matches1;
				if(matches1)
					System.out.println(line);
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			System.err.println("File "+args[0]+" doesn't exist!");
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new FilterTags(args);
	}

}
