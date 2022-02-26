import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Scanner;

import javax.imageio.ImageIO;


public class HiCStructure {
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 2)
		{
			System.err.println("Usaged: HiCStructure map1.txt map2.txt");
			System.err.println("Will output an image of the maps normalized by 5%-95% quntiles, a histogram of the interaction distances and a differential histogram, and the CDF of the distribution for various % of the map");
			System.exit(1);
		}
		int w = 0;
		int h = 0;
		double max1= 0;
		double min1= Double.POSITIVE_INFINITY;
//		double trueMax = 0;
		try{

			Scanner s1 = new Scanner(new File(args[0]));
			LinkedList<Double> l1 = new LinkedList<Double>();
			w=s1.nextLine().split("\t").length-2;
			while(s1.hasNextLine())
			{
				String[] line1 = s1.nextLine().split("\t");
				for(int i = 0; i < line1.length; i++)
				{
					if(i > 1)
					{
						double a= Double.parseDouble(line1[i]);
						l1.add(a);
					}
				}
				h++;
			}
			s1.close();
			Collections.sort(l1);
			min1 = l1.get(l1.size()/20);
			max1 = l1.get((int) (l1.size()*0.95));
			if(max1==min1) max1=min1+1;
			System.err.println("Min and max were: "+min1+"-"+max1);
			l1 = null;
			
			
		}catch(Exception e){e.printStackTrace();}


		try{
			BufferedImage img1 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			BufferedImage hist1 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			double histMax = 0;
			double sumHist = 0;
			double[] hist = new double[w];
			Scanner s1 = new Scanner(new File(args[0]));
			int j = 0;
			s1.nextLine();
			while (s1.hasNextLine() )
			{
				String[] line1 = s1.nextLine().split("\t");
				for(int i = 2; i < line1.length; i++)
				{
					double d = Double.parseDouble(line1[i]);
					img1.setRGB(i-2, j, getRGB(d,min1,max1));
					if(d > 0)
					{
						hist[Math.abs(i-2-j)]+=d;
						sumHist+=d;
						if(hist[Math.abs(i-2-j)] > histMax) histMax = hist[Math.abs(i-2-j)];
					}
				}
				j++;
				
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}
			System.err.println();
			Graphics2D g2 = (Graphics2D) hist1.getGraphics();
			g2.setColor(Color.white);
			g2.fillRect(0, 0, hist1.getWidth(), hist1.getHeight());
			g2.setColor(Color.red);
			for(int i = 1; i < 20; i++)
				g2.drawLine((i*w)/20, 0, (i*w)/20, h);
			g2.setColor(Color.black);
			//histMax/=2;
			for(int i = 0; i < hist.length; i++)
			{
				g2.drawLine(i,(int) (h-Math.round(h*(hist[i]/histMax))),i,h);
			}
			System.out.print("Index\t"+args[0]+"\t"+args[1]);
			
			//System.out.println();
			s1.close();
			ImageIO.write(img1, "png", new File(args[0]+".png"));
			System.err.println("Wrote image for "+args[0]);
			ImageIO.write(hist1, "png", new File(args[0]+"_hist.png"));
			System.err.println("Wrote histogram for "+args[0]);
			//now the idff
			w = 0;
			h = 0;
			max1= 0;
			min1= Double.POSITIVE_INFINITY;
			s1 = new Scanner(new File(args[1]));
			LinkedList<Double> l1 = new LinkedList<Double>();
			w=s1.nextLine().split("\t").length-2;
			while(s1.hasNextLine())
			{
				String[] line1 = s1.nextLine().split("\t");
				for(int i = 0; i < line1.length; i++)
				{
					if(i > 1)
					{
						double a= Double.parseDouble(line1[i]);
						l1.add(a);
					}
				}
				h++;
			}
			s1.close();
			Collections.sort(l1);
			min1 = l1.get(l1.size()/20);
			max1 = l1.get((int) (l1.size()*0.95));
			if(max1==min1) max1=min1+1;
			System.err.println("Min and max were: "+min1+"-"+max1);
			l1 = null;			
			img1 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			hist1 = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			double histMax2 = 0;
			double sumHist2 = 0;
			double[] hist2 = new double[w];
			s1 = new Scanner(new File(args[1]));
			j = 0;
			s1.nextLine();
			while (s1.hasNextLine() )
			{
				String[] line1 = s1.nextLine().split("\t");
				for(int i = 2; i < line1.length; i++)
				{
					double d = Double.parseDouble(line1[i]);
					img1.setRGB(i-2, j, getRGB(d,min1,max1));
					if(d > 0)
					{
						hist2[Math.abs(i-2-j)]+=d;
						sumHist2+=d;
						if(hist2[Math.abs(i-2-j)] > histMax2) histMax2 = hist2[Math.abs(i-2-j)];
					}


				}
				j++;
				
				if((j*100)/line1.length > ((j-1)*100)/line1.length)
					System.err.print(".");
			}
			System.err.println();
			g2 = (Graphics2D) hist1.getGraphics();
			g2.setColor(Color.white);
			g2.fillRect(0, 0, hist1.getWidth(), hist1.getHeight());
			g2.setColor(Color.red);
			for(int i = 1; i < 20; i++)
				g2.drawLine((i*w)/20, 0, (i*w)/20, h);
			g2.setColor(Color.black);

			s1.close();
			ImageIO.write(img1, "png", new File(args[1]+".png"));
			System.err.println("Wrote image for "+args[1]);
			ImageIO.write(hist1, "png", new File(args[1]+"_hist.png"));
			System.err.println("Wrote histogram for "+args[1]);
			System.out.println();
			for(int i = 0; i < hist.length; i++)
			{
				System.out.println((i+1)+"\t"+(hist[i]/sumHist)+"\t"+(hist2[i]/sumHist2));
			}
			//now lets output the differential histogram
			BufferedImage diffHist = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
			double maxDiffHist = 0;
			double minDiffHist = 999;
			for(int i = 0; i < Math.min(hist.length,hist2.length); i++)
			{
				hist[i]=(hist[i]/sumHist)-(hist2[i]/sumHist2);
				if(hist[i] > maxDiffHist)
					maxDiffHist=hist[i];
				if(hist[i] < minDiffHist)
					minDiffHist = hist[i];
			}
			g2 = (Graphics2D) diffHist.getGraphics();
			g2.setColor(Color.white);
			g2.fillRect(0, 0, diffHist.getWidth(), diffHist.getHeight());
			g2.setColor(Color.red);
			for(int i = 1; i < 20; i++)
				g2.drawLine((i*w)/20, 0, (i*w)/20, h);
			g2.drawLine(0,h/2,w,h/2);
			
			for(int i = 0; i < hist.length; i++)
			{
				if(hist[i] > 0)
				{
					g2.setColor(Color.black);
					g2.drawLine(i,(int) (h/2-Math.round((h/2)*(hist[i]/maxDiffHist))),i,h/2);
				}
				else
				{
					g2.setColor(Color.blue);
					g2.drawLine(i,(int) (h/2+Math.round((h/2)*(hist[i]/minDiffHist))),i,h/2);
				}
			}
			g2.setColor(Color.red);

			ImageIO.write(diffHist, "png", new File(args[0]+"-"+args[1]+"_hist.png"));
			System.err.println("Wrote histogram for "+args[0]+"-"+args[1]);
			
		}catch(Exception e){e.printStackTrace();}
		
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
	
}
