import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class JasparToHomer {

	private double getScore(int[] seq, double[][] pwm)
	{
		double score = 0;
		for(int i = 0; i < seq.length; i++)
		{
			score+=Math.log(pwm[seq[i]][i]/0.27);//using 0.27 instead of 0.25 to lower the threshold a bit. 
		}
		return score;
	}
	
	public JasparToHomer(String file) {
		try {
			Scanner s = new Scanner(new File(file));
			String header = s.nextLine();
			header=header+" Jaspar2020";
			String[] A=s.nextLine().split("\\s+");
			String[] C=s.nextLine().split("\\s+");
			String[] G=s.nextLine().split("\\s+");
			String[] T=s.nextLine().split("\\s+");
			double[][] counts = new double[4][A.length-3];
			int[] max = new int[counts[0].length];
			
			for(int i = 0; i < counts[0].length; i++)
			{
				double f= 0.0;
				double sum = 0;
				counts[0][i]=Double.parseDouble(A[2+i])+0.1;
				if(counts[0][i] > f)
				{
					f= counts[0][i];
					max[i]=0;
				}
				counts[1][i]=Double.parseDouble(C[2+i])+0.1;
				if(counts[1][i] > f)
				{
					f= counts[1][i];
					max[i]=1;
				}
				counts[2][i]=Double.parseDouble(G[2+i])+0.1;
				if(counts[2][i] > f)
				{
					f= counts[2][i];
					max[i]=2;
				}
				counts[3][i]=Double.parseDouble(T[2+i])+0.1;
				if(counts[3][i] > f)
				{
					f= counts[3][i];
					max[i]=3;
				}
				sum+=counts[0][i]+counts[1][i]+counts[2][i]+counts[3][i];
				counts[0][i]/=sum;
				counts[1][i]/=sum;
				counts[2][i]/=sum;
				counts[3][i]/=sum;
			}
			double threshold= 0;
			int tcount=0;
			for(int i = 0; i < counts[0].length; i++)
			{
				for(int j = 0; j < 4; j++)
				{
					int[] mut= max.clone();
					if(mut[i] != j)
					{
						mut[i] = j;
						threshold+=getScore(mut,counts);
						tcount++;
					}
				}
			}
			threshold/=tcount;
			System.out.println(header+"\t"+String.format("%3.4f",threshold));
			for(int i = 0; i < counts[0].length; i++)
			{
				System.out.println(String.format("%3.4f\t%3.4f\t%3.4f\t%3.4f",counts[0][i],counts[1][i],counts[2][i],counts[3][i]));
			}
			s.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args) {
		 		new JasparToHomer(args[0]);
	}

}
		
