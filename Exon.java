/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
/**
 * Defines a span of overlapping reads.
 * @author Maxim Shokhirev (C) 2015-2018
 *
 */
public class Exon implements Comparable<Exon>, Filterable
{
	int start, end;
	float count = 0;
	boolean plus;
	
	public Exon(int start, int end,boolean plus, float count)
	{
		this.start = start; this.end = end; this.plus = plus;this.count = count;
	}

	public double getDensity()
	{
		return (1.0*count)/Math.max(1,(end-start));
	}
	
	public int getScore()
	{
		return (int) Math.max(1,Math.min(1000,getDensity()*1000));
	}

	@Override
	public int compareTo(Exon o) {
		if(plus)
		{
			if(start < o.start)
			{
				if(end < o.start)
					return -2;   //start is less and not overlapping
				else
					return -1;   //start is less but still overlapping
			}
			else if(start > o.start)
			{
				if(start > o.end)
					return 2;  //start is greater and not overlapping
				else return 1;  //start is greater but still overlapping
			}
			else  //start the same
			{
				if(end < o.end)  //this one is shorter so it should go first
					return -1;
				else if(end > o.end)return 1;  
				else return 0;
			}
		}
		else
		{
			if(end < o.end)
			{
				if(end >= o.start)
					return 2;   //no overlap
				else
					return 1;   //overlap
			}
			else if(end > o.end) 
			{
				if(start > o.end)
					return -2;  //no overlap
				else return -1;  //overlap
			}
			else  //end the same
			{
				if(start > o.start)  //this one is shorter so it should go first
					return -1;
				else if(start < o.start)return 1;  
				else return 0;  //exactly the same exon
			}
		}
	}

	@Override
	public boolean filter() {
		if(Math.abs(end-start) >= 5 && count > 1)
			return false;
		return true;
	}
	@Override
	public String toString()
	{
		if(plus)
			return start+"-"+end+" + "+count;
		else return start+"-"+end+" - "+count;
	}
	
	public boolean overlaps(int start2, int end2, int ext)
	{
		double x = Math.min(end2, end);
		double y = Math.max(start2,start);
		if(x-y < -ext) return false;
		else return true;
	}
	
	
	public boolean overlaps(int start2, int end2)
	{
		double x = Math.min(end2, end);
		double y = Math.max(start2,start);
		if(x-y < 0) return false;
		else return true;
	}

	public double computeOverlapSize(Exon e) {
		double x = Math.min(e.end, end);
		double y = Math.max(e.start,start);
		if(x-y < 0) return 0;
		else return x-y;
	}
	
}
