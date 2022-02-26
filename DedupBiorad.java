import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Scanner;

public class DedupBiorad {

	HashSet<String> seen = new HashSet<String>();
	
	public DedupBiorad(String filename) {
		try {
			Scanner s = new Scanner(new File(filename));
			int count = 0;
			int total = 0;
			while(s.hasNextLine()){
			
				String line = s.nextLine();
				if(line.startsWith("@"))
				{
					System.out.println(line);
					continue;
				}
				String[] split = line.split("\t");
				String UMIplusPos = split[0].split(":")[0]+":"+split[2]+":"+split[3];
				if(!seen.contains(UMIplusPos))
				{
					seen.add(UMIplusPos);
					System.out.println(line);
					count++;
				}
				total++;
			}
			s.close();
			System.err.println("Kept "+count+"/"+total+" aligned reads from "+filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		if(args.length > 0)
		{
			new DedupBiorad(args[0]);
		}
		else System.err.println("Usage: DedupBiorad BioradWithUMIs.sam > Dedupped.sam");

	}

}
