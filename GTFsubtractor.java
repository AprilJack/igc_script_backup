import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


public class GTFsubtractor {

	public GTFsubtractor(String[] args) {
		if (args.length < 2)
		{
			System.err.println("Will subtract regions from a GTF based on another GTF. Only the last exon will be kept.");
			System.err.println("\nUsage: GTFsubtractor file1.gtf.gz file2.gtf.gz > file1-file2.gtf");
			System.exit(1);
		}
		GTF gtf1 = new GTF(args[0]);
		GTF gtf2 = new GTF(args[1]);
		int count =0;
		System.err.println("Outputing non-overlapping 3'UTR");
		for(String chrStr : gtf1.byChr.keySet())
		{
			boolean plus = chrStr.startsWith("+");
			for(Annotation a: gtf1.byChr.get(chrStr))
			{
				Region last = null;
				for(Region r: a.positions)
				{
					if(last == null || (plus && last.end < r.end)|| (!plus && last.start > r.start))
						last = r;
				}
				if(gtf2.byChr.get(chrStr) != null)
				{
					for(Annotation b: gtf2.byChr.get(chrStr)){
						if(b.overlaps(a) > 0)
						{
							//overlaps
							//lets just output the last exon of a trimmed down by the last exon of b...
							for(Region r: b.positions)
							{
								if(last.overlaps(r.start, r.end)>0 )
								{
									if(plus)
									{
										if(r.end >= last.start && r.end < last.end)
										{
											last.start = r.end;
										}
										else if(r.start > last.start && r.start <= last.end)
										{
											last.end = r.start;
										}
										else  //it completely covers
										{
											last.start=last.end;
											break;
										}
									}
									else
									{
										if(r.start > last.start && r.start <= last.end)
										{
											last.end = r.start;
										}
										else if(r.end >= last.start && r.end < last.end)
										{
											last.start = r.end;
										}
										else  //it completely covers
										{
											last.start=last.end;
											break;
										}
									}
								}
							}
							if(last.start == last.end)
								break;
						}
					}
				}
				//lets output the last region that is left
				if(last.start < last.end)
				{
					System.out.println(a.chr.substring(1)+"\tGTFSubtractor\ttranscript\t"+last.start+"\t"+last.end+"\t1000\t"+chrStr.charAt(0)+"\t.\t"+a.annotations);
					System.out.println(a.chr.substring(1)+"\tGTFSubtractor\texon\t"+last.start+"\t"+last.end+"\t1000\t"+chrStr.charAt(0)+"\t.\t"+a.annotations);
				}
				count++;
				if(count%1000==0) System.err.print(".");
			}
		}
		System.err.println();
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new GTFsubtractor(args);
	}

}
