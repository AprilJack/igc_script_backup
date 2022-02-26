import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class pos2start {

	public pos2start(String[] args) {
		if(args.length == 0)
		{
			System.err.println("Takes a HOMER position file and outputs a modified version centered on the TSS based on the strand.\n"
					+ "Usage: pos2start posFile > startPosFile");
		}
		else
		{
			try {
				Scanner s = new Scanner(new File(args[0]));
				while(s.hasNextLine())
				{
					String line = s.nextLine();
					String[] split = line.split("\t");
					if(split[4].compareTo("+")==0)
					{
						String out = split[0];
						for(int i = 1; i < split.length; i++)
						{
							if(i == 3)
								out+="\t"+split[2];
							else out+="\t"+split[i];
						}
						System.out.println(out);
					}
					else
					{
						String out = split[0];
						for(int i = 1; i < split.length; i++)
						{
							if(i == 2)
								out+="\t"+split[3];
							else out+="\t"+split[i];
						}
						System.out.println(out);
					}
					
				}
				s.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}

	public static void main(String[] args) {
		new pos2start(args);
	}

}
