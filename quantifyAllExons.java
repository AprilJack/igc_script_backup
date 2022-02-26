import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class quantifyAllExons {

	
	public quantifyAllExons(String[] args) {
		try {
			Scanner s = new Scanner(new File(args[0]));

			String current = "";
			
			double totalLength =0;
			int count = 0;
			String line =s.nextLine();
			String[] header = line.split("\t");
			double[] totalExp = new double[header.length-19];
			do
			{
				String[] split = line.split("\t");
				String id=split[0].split("-")[0];
				if(id.compareTo(current)!= 0 )
				{
					if(current.compareTo("")!= 0)
					{
						String out = current+"\t"+totalLength;
						for(int i = 0; i < totalExp.length; i++)
						{
							out+=String.format("\t%3.4f",1000.0*totalExp[i]/totalLength);
						}
						System.out.println(out);
						count++;
						if(count % 1000 == 0)System.err.print(".");
					}
					current=id; totalExp = new double[header.length-19];  totalLength= 0;
				}
				try{
					totalLength +=Math.abs(Double.parseDouble(split[3])-Double.parseDouble(split[2]))+1;
					for(int i = 19; i < split.length; i++)
						totalExp[i-19]+=Double.parseDouble(split[i]);
				}catch(Exception e)
				{
					System.err.println("Error on input: "+split.toString());
				}
				line = s.nextLine();
			}while(s.hasNextLine());
			System.err.println(" Finished!");
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new quantifyAllExons(args);
	}

}
