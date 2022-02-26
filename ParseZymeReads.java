
public class ParseZymeReads {

	private char[] template;
	private int[] map;
	private int[][] counts;
	private int end;
	private int length = -1;
	
	public ParseZymeReads(String[] args) {
		SuperScanner ss = new SuperScanner(args[0]);
		template = ss.getLine().split("  >  ")[0].toCharArray();
		ss.getLine(); //skip the alignment line
		end=Integer.parseInt(args[1]);

		if(args.length > 2)
			length=Integer.parseInt(args[2])-1;
		//0=match 1=mismatch A, 2=mismatch C, 3=mismatch G, 4=mismatch T, 5=deletion 6=insertion 7=insertion1 8=insertion2 ...26=insertion20 
		counts=new int[end][27];
		map = new int[template.length];
		int index = 0;
		map[0]=0;
		for(int i = 1; i < template.length; i++)
		{
			if(template[i]=='-')
				map[i]=index;
			else
			{
				index++;
				map[i]=index;
			}
		}
		int count =0;
		while(ss.hasMore())
		{
			String r = ss.getLine().split("  >  ")[0];
			processRead(r);
			count++;
			if(count%100000==0){
				System.err.println(count);
			}
		}
		//now let's output the results
		String header = "Position";
		if(length>-1)
			header=length+"Position";
		System.out.println(header+"\tMatches\tMismatch A\tMismatch C\tMismatch G\tMismatch T\tDeletions\tInsertions\tIns1\tIns2\tIns3\tIns4\tIns5\tIns6\tIns7\tIns8\tIns9\tIns10\tIns11\tIns12\tIns13\tIns14\tIns15\tIns16\tIns17\tIns18\tIns19\tIns20");
		for(int i = 0; i < counts.length; i++)
		{
			System.out.print(i+1);
			for (int j = 1; j < counts[i].length; j++)
			{
				System.out.print("\t"+counts[i][j-1]);
			}
			System.out.println();
		}
	}
	
	private void processRead(String read)
	{
		read=formatLeadingWhiteSpace(read);
		char[] c = read.toCharArray();
		if(length > -1)
		{
			int readLength= end;
			if(read.indexOf(' ') < map.length && read.indexOf(' ')>-1)
				readLength=map[read.indexOf(' ')-1];
			if(readLength != length)
				return;
		}
		int insertSize=0;
		for(int i = 0 ; i < Math.min(c.length,template.length); i++)
		{
			char t = template[i];
			 char r = c[i];
			if(r == ' '){
				break; // the read ended so let's stop counting.
			}
			if(c[i]!='-') 
			{
				if(t == r){
					counts[map[i]][0]++; //match
					if(insertSize>0)
					{
						counts[map[i-1]][6]++;
						counts[map[i-1]][insertSize+6]++;
						insertSize=0;
					}
				}
				else
				{
					if(t=='-') //insertion
					{
						insertSize++;
					}
					else
					{
						if(r=='A')
							counts[map[i]][1]++;//mismatch A
						else if(r=='C')
							counts[map[i]][2]++;//mismatch C
						else if(r=='G')
							counts[map[i]][3]++;//mismatch G
						else if(r=='T')
							counts[map[i]][4]++;//mismatch T
						if(insertSize>0)
						{
							counts[map[i-1]][6]++;
							counts[map[i-1]][insertSize+6]++;
							insertSize=0;
						}
					}
				}	
			}
			else //a dash in the read
			{
				if(t != r)//deletion
				{
					counts[map[i]][5]++;
					if(insertSize>0)  //an insert that ends with a deletion
					{
						counts[map[i-1]][6]++;
						counts[map[i-1]][insertSize+6]++;
						insertSize=0;
					}
				}
			}
		}
	}

	private String formatLeadingWhiteSpace(String read) {
		// TODO Auto-generated method stub
		char[] chars = read.toCharArray();
		for(int i = 0; i < chars.length; i++)
		{
			if(chars[i]==' ')
				chars[i]='-';
			else break;
		}
		return new String(chars);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length > 1)
			new ParseZymeReads(args);
		else
			System.err.println("Usage: ParseZymeReads breseqBamToAln TemplateLength");
	}

}
