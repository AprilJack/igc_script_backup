import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.zip.GZIPOutputStream;


public class gtfFixer {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SuperScanner s = new SuperScanner(args[0]);
		boolean needsFixin = false;
		while(s.hasMore())
		{
			String line = s.getLine();
			String[] split = line.split("\t");
			if(split[2].compareTo("transcript")==0 && split[3].contains("-"))
			{
				needsFixin = true;
				break;
			}
		}
		s.close();

		if(needsFixin)
		{
			try {
				PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(new File(args[0]+"_tmp"))));
				s = new SuperScanner(args[0]);
			
				while(s.hasMore())
				{
					String line = s.getLine();
					String[] split = line.split("\t");
					if(split[2].compareTo("transcript")==0 && split[3].contains("-"))
					{
						pw.println(split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3].split("-")[0]+"\t"+split[3].split("-")[1]+"\t"+split[4]+"\t"+split[5]+"\t"+split[6]+"\t"+split[7]);
					}
					else pw.println(line);
				}
				pw.close();
				s.close();
				File f = new File(args[0]+"_tmp");
				File old =new File(args[0]);
				old.delete();
				f.renameTo(old);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.err.println("Done fixing "+args[0]);
		}
		else
			System.err.println(args[0]+" doesn't need fixing.");
	}

}
