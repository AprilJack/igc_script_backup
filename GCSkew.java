import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;


public class GCSkew {

	String seq = "";
	int wsize = 100;
	LinkedList<Byte> window = new LinkedList<Byte>();
	int C = 0;
	int g = 0;
	int a = 0;
	int t = 0;
	public GCSkew(String[] args) {
		try {
			if(args.length < 1)
			{
				System.err.println("Usage: GCSkew sequenceFile.txt\nOutputs 4 columns: index, GCSkew, ATSkew, base to stdout for plotting. Non-AaTtGgCc are ignored");
				System.exit(1);
			}
			Scanner s = new Scanner(new File(args[0]));
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				seq+=line;
			}
			s.close();
			int count = 1;
			char[] sequence = seq.toCharArray();
			System.out.println("Index\tGCSkew\tATSkew\tBase");
			for(char c: sequence )
			{
				if(c == 'A' || c == 'a')
				{
					a++;
					window.add((byte) 0);
				}
				else if(c == 'C' || c == 'c')
				{
					C++;
					window.add((byte) 1);
				}
				else if(c == 'G' || c == 'g')
				{
					g++;
					window.add((byte) 2);
				}
				else if(c == 'T' || c == 't')
				{
					t++;
					window.add((byte) 3);
				}
				if(wsize < window.size())
				{
					byte first = window.removeFirst();
					if(first == 0)
						a--;
					else if(first == 1)
						C--;
					else if(first == 2)
						g--;
					else t--;
				}
				double gcskew = 0;
				if(C+g > 0)
					gcskew = (C-g)/(1.0*C+g);
				double atskew = 0;
				if(a+t > 0)
					atskew = (a-t)/(1.0*a+t);
				System.out.println(String.format("%d\t%3.3f\t%3.3f\t%s",count,gcskew,atskew,c));
				count++;
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		new GCSkew(args);

	}

}
