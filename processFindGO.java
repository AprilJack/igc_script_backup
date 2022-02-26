import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;


public class processFindGO {
	HashMap<String,ArrayList<Term>> map = new HashMap<String,ArrayList<Term>>();
	String[] databases = {"kegg","wikipathways","biological_process","pfam","cosmic","msigdb","gwas","reactome","biocyc"};

	
	private class Term implements Comparable<Term>
	{
		double q;
		int count = 0;
		String genes;
		String name;
		String source;
		
		public Term(String name, String genes, double q, String source)
		{
			this.q=q;this.name=name; this.source = source;
			if(genes.length() > 32000)
				genes=genes.substring(0,32000);
			this.genes =genes;
			sortGenes();
		}
		
		private void sortGenes()
		{
			LinkedList<String> temp = new LinkedList<String>();
			String[] split = genes.split(",");
			for(int i = 0; i < split.length; i++)
				temp.add(split[i]);
			Collections.sort(temp);
			genes = "";
			count = temp.size();
			for(String s: temp)
			{
				if(genes.length() == 0)
					genes = s;
				else
					genes+=","+s;
			}
		}

		@Override
		public int compareTo(Term o) {
			// TODO Auto-generated method stub
			if(o.q < q ) return 1;
			else if(o.q > q) return -1;
			else return 0;
		}
		
		@Override
		public String toString()
		{
			return name+"("+q+")";
		}
		
	}
	
			
	public processFindGO(String[] args) {
		if(args.length > 1)
		{
			databases = new String[args.length-1];
			for(int i = 1; i < args.length; i++)
				databases[i-1]=args[i];
		}
		File dir = new File(args[0]);
		processDir(dir);
		//lets output the map
		String header = "";
		for(String file: map.keySet())
		{
			Collections.sort(map.get(file));
		}
		for(String file: map.keySet())
		{
			String flabel = file;
			if(file.indexOf(File.separatorChar) != -1)
				flabel = file.substring(file.lastIndexOf(File.separatorChar)+1);
			if(header.length() == 0)
			{
				header+="Terms\t"+flabel+"\t#Genes\tGenes\tSource";
			}
			else
				header+="\tTerms\t"+flabel+"\t#Genes\tGenes\tSource";
		}
		System.out.println(header);

		int h = 0;

		while(true)
		{
			try{
				String s = "";
				for(String file: map.keySet())
				{
					Term t = map.get(file).get(h);
					if(s.length() == 0)
						s+=t.name+"\t"+t.q+"\t"+t.count+"\t"+t.genes+"\t"+t.source;
					else s+="\t"+t.name+"\t"+t.q+"\t"+t.count+"\t"+t.genes+"\t"+t.source;
				}
				System.out.println(s);
				h++;
			}catch(Exception e){break;}
		}
		
		
	}

	private void processDir(File dir)
	{

		for(File f: dir.listFiles())
		{
			if(f.isDirectory())
			{
				processDir(f);
				System.err.println("Processed "+f.getName());
			}
			else
			{
				String db = null;
				for(String database: databases)
				{
					if(f.getName().contains(database))
					{
						db =f.getName();
						break;
					}
				}
				if(db != null)
				{
					//lets try to process it!
					try {
						System.err.println("Found file "+f.getName());
						HashSet<String> lines = new HashSet<String>();
						Scanner s = new Scanner(f);
						if(s.hasNextLine())
						s.nextLine();  //header
						//LinkedList<Double> pvalues = new LinkedList<Double>();
						while(s.hasNextLine())
						{
							String sline =s.nextLine();
							String[] split = sline.split("\t");
							boolean found = false;
							for(String line: lines)
							{
								if(line.contains(split[1])||split[1].contains(line)) //checking if this term already found before (e.g. kegg and wikipathways have similar terms)
								{
									found = true;
									break;
								}
							}
							try{
								if(!found)
								{
									if(split.length > 10)
										lines.add(split[1]+"\t"+split[2]+"\t"+split[10]);
									else
										lines.add(split[1]+"\t"+split[2]+"\t ");
								}
							} catch(Exception e){System.err.println("Encountered poorly formatted GO term :"+sline);}
						}
						LinkedList<Term> terms = new LinkedList<Term>();
						for(String str: lines)
						{
							//found the category of interest
							double pvalue = Double.parseDouble(str.split("\t")[1]);
							//pvalues.add(pvalue);
							if(map.get(f.getParent()) == null)
								map.put(f.getParent(), new ArrayList<Term>());
							Term newTerm = new Term(str.split("\t")[0],str.split("\t")[2],pvalue,f.getName());
							//map.get(f.getParent()).add(newTerm);
							terms.add(newTerm);
						}
						int count = terms.size();
						double lowestQThusFar=1.0;
						//ok now lets adjust the p-values using the Storey (2002) approach:
						//A direct approach to false discovery rates. John D. Storey. J. R. Statist. Soc
						//Iterator<Double> piter = pvalues.iterator();  //these will be in order of smallest to largest
						//lets assume that pi0 is just 1.0 (e.g. we suspect that the vast majority of tests will come from the null)
						Collections.sort(terms);
						Collections.reverse(terms);
						for(Term t: terms)
						{
							t.q = Math.min(lowestQThusFar,(t.q*terms.size())/count);
							lowestQThusFar=Math.min(t.q,lowestQThusFar);
							count--;
						}
						if(map.get(f.getParent()) == null)
							map.put(f.getParent(), new ArrayList<Term>());
						map.get(f.getParent()).addAll(terms);
						s.close();


					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				}
		
			}
		}
	}

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length < 1)
		{
			System.err.println("Usaged: processFindGO Dir_with_comparison_dirs [database1] [database2] [database3], ...");
			System.err.println("\nExample: processFindGO . kegg msigdb biological_process cosmic pfam wikipathways reactome biocyc\n defaults to these.");
			System.err.println("Prints a tab file to STDOUT containing enriched terms, FDRs, and genes");
		}
		else
			new processFindGO(args);
	}

}
