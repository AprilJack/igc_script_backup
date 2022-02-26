import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.zip.GZIPInputStream;


public class FastScanner {

	protected byte[] buffer = new byte[1024];
	protected int pos = 0;
	protected int max = 0;
	protected boolean open = true;
	protected BufferedInputStream b = null;
	
	public boolean hasMore()
	{
		return (open && max > 0);
	}
	
	public FastScanner(String file)
	{
		try {
			if(!file.endsWith(".gz"))
			{
				b = new BufferedInputStream(new FileInputStream(new File(file)),8096);
				max = b.read(buffer);
			}
			else
			{
				b = new BufferedInputStream(new GZIPInputStream(new FileInputStream(new File(file))),8096);
				max = b.read(buffer);
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public String getLine()
	{
		if(!open) return null;
		StringBuilder sb = new StringBuilder(1024);
		while(true)
		{
			while(pos < max)
			{
				char c = (char) buffer[pos++];
				if (c == '\n')
				{
					String line = sb.toString();
					sb = new StringBuilder(1024);
					return line;
				}
				else if(c != '\r')
				{
					sb.append(c);
				}
			}
			
			try {
				max = b.read(buffer);
				pos = 0;  //rewind
				if(max == -1)
				{
					open = false;
					close();
					return sb.toString();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public void close()
	{
		try {
			b.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
