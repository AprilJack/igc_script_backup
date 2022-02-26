import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


public class FastaVsFasta {

	HashMap<String,HashSet<String>> map = new HashMap<String,HashSet<String>>();
	ArrayList<String> files = new ArrayList<String>();
	public FastaVsFasta(String[] args) {
		if(args.length < 1)
			System.err.println("Usage: a.fasta b.fasta [c.fasta] ... [n.fasta]\nProduces stats on all pairwise comparisons of sequences.");
		for(String s: args)
		{
			map.put(s, new HashSet<String>());
			SuperScanner ss = new SuperScanner(s);
			files.add(s);
			while(ss.hasMore())
			{
				ss.getLine();
				map.get(s).add(ss.getLine());
			}
			ss.close();
		}
		for(int i = 0; i < files.size(); i++)
		{
			for(int j = i+1; j < files.size(); j++)
			{
				System.out.println("Uniq to "+files.get(i)+"\tUniq to "+files.get(j)+"\tIn both\tTotal");
				int[] counts = compare(map.get(files.get(i)), map.get(files.get(j)));
				System.out.println(String.format("%d\t%d\t%d\t%d",counts[0],counts[1],counts[2],counts[0]+counts[1]+counts[2]));
			}
		}
	}
	

	private int[] compare(HashSet<String> a, HashSet<String> b) {
		int[] result = new int [3];
		HashSet<String> both = new HashSet<String>();
		both.addAll(a); both.addAll(b);
		for(String s: both)
		{
			if(a.contains(s))
			{
				if(b.contains(s))
					result[2]++;
				else
					result[0]++;
			}
			else
				result[1]++;
		}
		return result;
	}


	public static void main(String[] args) {
		new FastaVsFasta(args);

	}

}
