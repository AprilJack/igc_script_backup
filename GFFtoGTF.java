import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Scanner;


public class GFFtoGTF {

	public GFFtoGTF(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: GFFtoGTF GFF3file feature1 [feature2] ... > output.gtf\n"
					+ "ex: GFFtoGTF annotation.gff exons transposable_elements tRNA > output.gtf\n"
			);
			System.exit(1);
		}
		HashSet<String> features = new HashSet<String>();
		for(int i = 1; i < args.length; i++)
			features.add(args[i]);
		try {
			Scanner s = new Scanner(new File(args[0]));
			while(s.hasNextLine())
			{
				String line = s.nextLine();
				String[] split = line.split("\t");
				if(split.length > 1)
				{
					if(features.contains(split[2]))
					{
						//feature of interest... lets output in GTF format!
						String[] info = split[8].split(";");
						String gene_id = "";
						String description = "";
						for(String str: info)
						{
							if(str.startsWith("ID=")||str.startsWith("Parent="))
								gene_id = str.split("=")[1].split(":")[0].trim();
							else
							{
								//description+= str.split("=")[0].trim()+"=\""+str.split("=")[1].trim()+"\"; ";
								description+= str;
							}
						}
						System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t.\t"+split[6]+"\t0.0\tgene_id \""+gene_id+"\"; transcript_id \""+gene_id+"\"; "+description);
					}
				}
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new GFFtoGTF(args);
	}

}
