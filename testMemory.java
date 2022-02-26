public class testMemory
{
	public static void main (String[] args)
	{
		System.err.println("Starting memory usage: "+((Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024)));
		byte[] test = new byte[100000000*4];
		System.err.println("Ending memory usage: "+((Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024)));
		test = null;
		System.gc();
		byte[][] test22 = new byte[100000000][4];
		System.err.println("Ending memory usage: "+((Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024)));
	}


}
