import java.util.LinkedList;


public class nmerGenerator {
	
	private void generate2(int n, LinkedList<String> nmers,String nmer)
	{
		if(n == 0)
		{
			nmers.add(nmer);
		}
		else
		{
			generate2(n-1,nmers,nmer+"A");
			generate2(n-1,nmers,nmer+"C");
			generate2(n-1,nmers,nmer+"G");
			generate2(n-1,nmers,nmer+"T");
		}
	}

	public nmerGenerator(int n) {
		LinkedList<String> nmers = new LinkedList<String>();
		generate(n,nmers);
		for(String nmer: nmers)
		{
			System.out.println(">"+nmer);
			System.out.println(nmer);
		}
	}

	private void generate(int n, LinkedList<String> nmers) {
		generate2(n-1,nmers,"A");
		generate2(n-1,nmers,"C");
		generate2(n-1,nmers,"G");
		generate2(n-1,nmers,"T");		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new nmerGenerator(Integer.parseInt(args[0]));
	}

}
