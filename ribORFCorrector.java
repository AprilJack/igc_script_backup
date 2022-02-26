import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;


public class ribORFCorrector {

	public ribORFCorrector(String[] args) {
		if(args.length < 2)
		{
			System.err.println("ribORFCorrector length read.dist.sample.length >> offset.correction.parameters");
			System.err.println("Assumes your read.dist.sample.length file is from -30 ");
			System.exit(1);
		}
		File f = new File(args[1]);
		if(f.exists())
		{
			try {
				Scanner s = new Scanner(f);
				String[] start = s.nextLine().split("\t");
				double[] mean = new double[3];
				double[][] val = new double[(start.length-1)/3+1][3];
				
				for(int i = 1; i < start.length; i++)
				{
					mean[(i-1)%3]+=Double.parseDouble(start[i]);
					val[(i-1)/3][(i-1)%3]=Double.parseDouble(start[i]);
				}
				String[] end = s.nextLine().split("\t");
				for(int i = 1; i < end.length; i++)
				{
					mean[i%3]+=Double.parseDouble(end[i]);
				}
				if(mean[0]*0.8>mean[1] && mean[0]*0.8>mean[2])
				{
					//good quality
					double m = mean[0]/((2*start.length-1)/3);
					for(int i =0; i < val.length; i++)
					{
						if(val[i][0]>m*1.5)
						{
							int shift = 30-(i*3)+3;
							System.out.println(args[0]+"\t"+shift);
							break;
						}
					}
				}
				else if(mean[1]*0.8>mean[0] && mean[1]*0.8>mean[2])
				{
					//good quality
					double m = mean[1]/((2*start.length-1)/3);
					for(int i =0; i < val.length; i++)
					{
						if(val[i][1]>m*1.5)
						{
							int shift = 30-(i*3)+2;
							System.out.println(args[0]+"\t"+shift);
							break;
						}
					}
				}
				else if(mean[2]*0.8>mean[0] && mean[2]*0.8>mean[1])
				{
					//good quality
					double m = mean[2]/((2*start.length-1)/3);
					for(int i =0; i < val.length; i++)
					{
						if(val[i][2]>m*1.5)
						{
							int shift = 30-(i*3)+1;
							System.out.println(args[0]+"\t"+shift);
							break;
						}
					}
				}
				s.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		else System.err.println("File "+args[1]+ " does not exist!");
	}
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ribORFCorrector(args);
	}

}
