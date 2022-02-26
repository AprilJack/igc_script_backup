
public class Region {
	int start, end;
	public Region(int start, int end)
	{
		this.start=start;this.end=end;
	}
	
	public int overlaps(int left, int right)
	{
		return Math.max(0,Math.min(right,end)-Math.max(left,start));
	}
	
	public int size()
	{
		return end-start;
	}
}
