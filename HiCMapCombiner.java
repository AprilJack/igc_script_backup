import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Scanner;

import javax.imageio.ImageIO;


public class HiCMapCombiner {
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 2)
		{
			System.err.println("Usaged: HiCMapCombiner map1.txt map2.txt [smoothing window=0] [threshold=25]");
			System.err.println("Will output binary comparison between map1 and map2 average withing the given window. Extreme images are colored by threshold");
			System.exit(1);
		}
		int w = 0;
		int h = 0;
		int win =0;
		int threshold = 25;
		double max1= 0;
		double max2 =0;
		double min1= Double.POSITIVE_INFINITY;
		double min2= Double.POSITIVE_INFINITY;
		double best_m = 0;
		if(args.length > 2)
		{
			win = Integer.parseInt(args[2]);
		}
		if(args.length > 3)
		{
			threshold = Math.min(255,Integer.parseInt(args[3]));
			
		}

		//lets try to sample the data
		double avg1 = 0;
		double avg2 = 0;


		try{
			Scanner s1 = new Scanner(new File(args[0]));
			Scanner s2 = new Scanner(new File(args[1]));
			LinkedList<Double> l1 = new LinkedList<Double>();
			LinkedList<Double> l2 = new LinkedList<Double>();
			w=s1.nextLine().split("\t").length-2; s2.nextLine(); //ignore header
			while(s1.hasNextLine() && s2.hasNextLine())
			{
				String[] line1 = s1.nextLine().split("\t");
				String[] line2 = s2.nextLine().split("\t");
				for(int i = 0; i < Math.min(line1.length,line2.length); i++)
				{
					if(i > 1)
					{
						double a= Double.parseDouble(line1[i]);
						double b= Double.parseDouble(line2[i]);
						avg1 += a;
						avg2 += b;
						l1.add(a);
						l2.add(b);
					}
				}
				h++;
			}
			s1.close();
			s2.close();
			Collections.sort(l1);
			Collections.sort(l2);
			min1 = l1.get(l1.size()/20);
			min2 = l2.get(l2.size()/20);
			max1 = l1.get((int) (l1.size()*0.95));
			max2 = l2.get((int) (l2.size()*0.95));
			if(max1==min1) max1=min1+1;
			if(max2==min2) max2=min2+1;
			System.err.println(min1+"-"+max1+" "+min2+"-"+max2);
			//System.err.println(l1.getFirst()+" "+l1.getLast()+" "+l2.getFirst()+" "+l2.getLast());
			l1 = null;
			l2 = null;
			
		}catch(Exception e){e.printStackTrace();}
		best_m = avg1/avg2;

		System.err.println(args[0]+" was "+best_m+" x "+args[1]);

		StringBuilder sb = new StringBuilder(1024*1024*128);
		double[][] combined = new double[w][h];

		try{
			BufferedImage img1 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			BufferedImage img2 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			BufferedImage img3 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			BufferedImage img3b = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			Scanner s1 = new Scanner(new File(args[0]));
			Scanner s2 = new Scanner(new File(args[1]));
			sb.append(s1.nextLine());  //assume headers are the same
			sb.append(System.lineSeparator());
			s2.nextLine();
			double diff = 0;
			int count = 0;
			int j = 0;
			while (s1.hasNextLine() && s2.hasNextLine())
			{
				String[] line1 = s1.nextLine().split("\t");
				String[] line2 = s2.nextLine().split("\t");
				for(int i = 0; i < line1.length; i++)
				{

					if (i == 0)sb.append(line1[i]);
					else if(i ==1)sb.append("\t"+line1[i]);
					else
					{
						double d1 = Double.parseDouble(line1[i]);
						double d2 = Double.parseDouble(line2[i]);
						if(d1 >=max1 && d2 <max2)
							combined[i-2][j]= -1;
						else if(d1 < max1 && d2 >= max2)
							combined[i-2][j]= 1;
						//double d3 = (d1-d2*best_m);
						//combined[i-2][j]=d3;
						sb.append(String.format("\t%3.3f",combined[i-2][j]));

						img1.setRGB(i-2, j, getRGB(d1,min1,max1));
						img2.setRGB(i-2, j, getRGB(d2,min2,max2));
					}
				}
				sb.append(System.lineSeparator());
				j++;
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}
			double low = min1-max2*best_m;
			double high = max1-min2*best_m;
			double cmax = 0;
			double cmin = 0;
			int totalMinus = 0;
			int totalPlus =0;
			for(int i = 0; i < combined.length; i++)
			{
				for(j = 0; j < combined[i].length; j++)
				{
					double avg = 0;
					int ct = 0;
					for(int x=-win; x<=win; x++)
					{
						for(int y=-win; y <=win; y++)
						{
							if(i+x >=0 && i+x < combined.length && j+y >=0 && j+y < combined[i].length)
							{
								double toAverage =combined[i+x][j+y];
								avg= avg+toAverage;
								ct++;
							}
						}
					}
					avg = avg/ct;
					double pavg = avg;
					if(avg < 0) pavg = -avg;
					diff += pavg;
					if(avg < cmin) cmin = avg;
					if(avg > cmax) cmax = avg;
					//System.err.println(diff+" "+avg+" "+ct+" "+count);
					count++;
					img3.setRGB(i, j, getRGB(avg,low,high));
					img3b.setRGB(i, j, getExtremeRGB(avg,low,high,threshold));
				}
			}
			String fraction_different = String.format("%3.2f",(100.0*(totalMinus+totalPlus))/(w*h));
			System.err.println("Done");
			System.err.println("Average % diff was :"+fraction_different+" (diff/count)="+(diff/count)+" cmin="+cmin+" cmax="+cmax);
			ImageIO.write(img1, "png", new File(args[0]+".png"));
			System.err.println("Wrote image for "+args[0]);
			ImageIO.write(img2, "png", new File(args[1]+".png"));
			System.err.println("Wrote image for "+args[1]);
			ImageIO.write(img3, "png", new File(combineFileNames(args[0],args[1])+".png"));
			ImageIO.write(img3b, "png", new File(combineFileNames(args[0],args[1])+"_"+fraction_different+"_extreme.png"));
			System.err.println("Wrote image combined subtracted map");
			System.out.print(sb.toString());
			s1.close(); s2.close();

		}catch(Exception e){e.printStackTrace();}
		
	}
	
	private static String combineFileNames(String f1, String f2) {
		String result = f1;
		if(f2.indexOf(File.pathSeparatorChar)!= -1)
		//if(f2.indexOf("/")!= -1)
		{
			result+="-"+f2.substring(f2.lastIndexOf(File.pathSeparatorChar)+1);
			//result+="-"+f2.substring(f2.lastIndexOf("/")+1);
		}
		else
			result+="-"+f2;
		return result;
	}

	private static int getRGB(double d,double min, double max)
	{
		if(d <= 0)
		{
			int c = (int) Math.min(255,255*(1.0-(d-min)/(-min)));

			return  (((255-c)&0x0ff)<<16)|((255&0x0ff)<<8)|((255-c)&0x0ff);
		}
		else
		{
			int c = (int) Math.min(255,(255*(d))/(max));
			return  (((255)&0x0ff)<<16)|(((255-c)&0x0ff)<<8)|((255-c)&0x0ff);
		}
	}
	private static int getExtremeRGB(double d,double min, double max, int threshold)
	{
		if(d <= 0)
		{
			int c = (int) Math.min(255,255*(1.0-(d-min)/(-min)));
			if(c < threshold)
				c=0;
			else
				c=255;
			return  (((255-c)&0x0ff)<<16)|((255&0x0ff)<<8)|((255-c)&0x0ff);
		}
		else
		{
			int c = (int) Math.min(255,(255*(d))/(max));
			if(c < threshold)
				c=0;
			else
				c=255;
			return  (((255)&0x0ff)<<16)|(((255-c)&0x0ff)<<8)|((255-c)&0x0ff);
		}
	}
	
}
