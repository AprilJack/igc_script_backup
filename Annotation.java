
import java.util.LinkedList;



public class Annotation {
	LinkedList<Region> positions = new LinkedList<Region>();
	boolean coding = false;  //does it have a CDS annotation?
	String chr = "";
	String id = "";
	String annotations;
	
	public Annotation(String chr, String id)
	{
		this.chr = chr; this.id = id;
	}
	
	public void addRegion(int start, int end)
	{
		positions.add(new Region(start,end));
	}
	
	public LinkedList<Region> getRegions()
	{
		return positions;
	}
	
	public void addAnnotation(String str)
	{
		annotations = str;
	}
	
	
	public int overlaps(Annotation a)
	{
		if(Math.min(a.getMax(), getMax())-Math.max(a.getMin(),getMin())> 0)
			return Math.min(overlapsAgainst(a), a.overlapsAgainst(this));
		else
			return 0;
	}
	
	/**
	 * Returns a value corresponding to the quality of the overlap
	 * @param a
	 * @return 0 = no overlap at all, 
	 * 1 = some exons are overlapped
	 * 2 = some exons are perfectly overlapped
	 * 3 = all exons are overlapped
	 * 4 = internal exons are perfectly matched
	 * 5 = all exons are perfectly matched (identical)
	 */
	private int overlapsAgainst(Annotation a)
	{	
		int overlapped = 0;
		int matched = 0;
		int internalMatched = 0;
		int rcount = 0;
		for(Region r1: positions)
		{
			int bestOverlap = 0;
			boolean startMatched = false;
			boolean endMatched = true;
			for(Region r2: a.positions)
			{
				int overlap = r1.overlaps(r2.start,r2.end);
				if(overlap > bestOverlap)
				{
					bestOverlap = overlap;
					if(r1.start == r2.start) startMatched =true;
					if(r1.end == r2.end) endMatched = true;
				}
			}
			if(bestOverlap > 0)
			{
				if(bestOverlap != r1.size())
					overlapped++;
				else{
					matched++;
				}
				if((rcount==0 && endMatched)||(rcount==positions.size()-1 && startMatched)||(startMatched&&endMatched))
					internalMatched++;
			}
			rcount++;
		}
		if(matched == positions.size())
			return 5;
		else if(internalMatched == positions.size())
			return 4;
		else if(overlapped == positions.size())
			return 3;
		else if(matched > 0)
			return 2;
		else if(overlapped > 0)
			return 1;
		else return 0;
	}
	
	public int getTranslatedLength()
	{
		int count = 0;
		for(Region r: positions)
		{
			count+=r.end-r.start;
		}
		return count;
	}
	
	public String toString()
	{
		String codeString = "Coding";
		if(!coding)codeString ="Non-coding/Unknown";
		return id+"\t"+codeString+"\t"+chr+"\t"+positions.getFirst().start+"\t"+positions.getLast().end+"\t"+annotations;
	}

	public int getMin() {
		return Math.min(positions.getFirst().start,positions.getLast().start);
	}
	
	public int getMax() {
		return Math.max(positions.getFirst().end,positions.getLast().end);
	}
}
