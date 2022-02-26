import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class gbk2gtf {
	
	public gbk2gtf(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: gbk2gtf file.gbk chr_name > file.gtf");
			System.exit(1);
		}
		try {
			Scanner s  = new Scanner(new File(args[0]));
			String line="";
			String sign = "+";
			String geneID = null;
			String note = null;
			String[] locations = null;
			String feature="";
			while(s.hasNextLine())
			{
				if(line.length() == 0)
					line =s.nextLine();

				if(line.startsWith("     ") && (line.startsWith("     CDS")||line.startsWith("     misc_RNA")||line.startsWith("     exon")||
						line.startsWith("     repeat_region")||line.startsWith("     enhancer")||line.startsWith("     splicing_signal")||line.startsWith("     CAAT_signal")
						||line.startsWith("     polyA_site"))||line.startsWith("     misc_marker")||line.startsWith("     TATA_signal"))
				{
					//first lets check if we want to output previous feature
					if(locations != null)
					{
						if(feature.compareTo("CDS")==0)
							feature = "exon";
						//lets output everything
						if(geneID== null)geneID=note;
						for(String range: locations)
						{
							String[] ranges = range.split("\\.\\.");
							System.out.println(String.format("%s\tgbk2gtf\t%s\t%s\t%s\t.\t%s\t0\tgene_id \"%s\"; transcript_id \"%s\"; notes \"%s\";",args[1],feature,ranges[0],ranges[1],sign,geneID,geneID,note));
						}
						note = null;
						locations = null;
						sign="+";
						geneID = null;
						feature = null;
					}
					line=line.replaceAll("join\\(", "");
					line=line.replaceAll("\\)","");
					line=line.replaceAll(">","");
					line=line.replaceAll("<","");
					if(line.contains("complement"))
					{
						line = line.replaceAll("complement\\(", "");
						sign = "-";
					}
					feature = line.split("\\s+")[1];
					locations= line.split(feature)[1].trim().split(",");
					line=s.nextLine();
				}
				else
				{
					line = s.nextLine();
					if(line.contains("="))
					{
						if(line.contains("product="))
						{
							geneID =line.split("=")[1].replaceAll("\"","");
						}
						else
						{
							String info = line.split("=")[1].replaceAll("\"","");
							if(note == null) note = info;
							else note+=","+info;
						}
					}
				}

			}
			if(locations != null)
			{
				if(feature.compareTo("CDS")==0)
					feature = "exon";
				//lets output everything
				if(geneID== null)geneID=note;
				for(String range: locations)
				{
					String[] ranges = range.split("\\.\\.");
					System.out.println(String.format("%s\tgbk2gtf\t%s\t%s\t%s\t.\t%s\t0\tgene_id \"%s\"; transcript_id \"%s\"; notes \"%s\";",args[1],feature,ranges[0],ranges[1],sign,geneID,geneID,note));
				}
				note = null;
				locations = null;
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args)
	{
		new gbk2gtf(args);
	}

}
