/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
import java.util.HashSet;
import java.util.LinkedList;

/**
 * Wrapper for a collection of Exons
 * @author Maxim Shokhirev (C) 2016
 *
 */
public class Transcript {
	private LinkedList<Exon> exons = new LinkedList<Exon>();
	private double TPM = 0;
	private double FPKM = 0;
	private double error = 0;
	
	public double getError() {
		return error;
	}

	public void setError(double error) {
		this.error = error;
	}

	public LinkedList<Exon> getExons()
	{
		return exons;
	}
	
	public Exon get(int i)
	{
		return exons.get(i);
	}
	
	public void add(Exon e)
	{
		exons.add(e);
	}
	
	public void setTPM(double TPM)
	{
		this.TPM = TPM;
	}
	
	public void setFPKM(double FPKM)
	{
		this.FPKM = FPKM;
	}

	public Transcript duplicate() {
		Transcript dup = new Transcript();
		dup.exons.addAll(exons);
		dup.TPM = TPM;
		dup.FPKM = FPKM;
		return dup;
	}
	/**
	 * the minimum coord
	 * @return
	 */
	public int getStart() {
		return Math.max(1,Math.min(exons.getFirst().start,exons.getLast().start));
	}
	/**
	 * The maximum coord
	 * @return
	 */
	public int getEnd() {
		return Math.max(exons.getFirst().end,exons.getLast().end);
	}
	
	public boolean isPlus()
	{
		if(exons.size()==0) return true;
		return exons.get(0).plus;
	}

	public double getAvgIntensity() {
		double avgIntensity = 0;
		for(Exon e: exons)
		{
			avgIntensity+=e.getDensity();
		}
		return (avgIntensity/=exons.size());
	}
	
	public int getScore()
	{
		return (int) Math.max(1,Math.min(1000,getAvgIntensity()*1000));
	}
	/**
	 * Returns true if this transcript cover at least fraction of t
	 * @param t
	 * @param fraction
	 * @return
	 */
	public boolean covers(Transcript t, double fraction)
	{
		//first lets check if they cover about the same length
		double size = t.getEnd()-t.getStart();   //size of the transcript that this one is being compared to
		double overlap = Math.min(getEnd(),t.getEnd())-Math.max(getStart(),t.getStart()); //overlap length
		if(overlap > fraction*size)  //this transcript overlaps at least fraction of t
			return true;  //they overlap on the genome and are about the same size
		else return false;
	}
	
	/**
	 * Returns true if there is any overlap between exonic bases
	 * @param t
	 * @return
	 */
	public boolean exonsCover(Transcript t)
	{
		if(covers(t,0))
		{
			for(Exon e: exons)
			{
				for(Exon e2: t.exons)
				{
					if(e.overlaps(e2.start-1, e2.end+1)) //to ensure no intron exons
						return true;
				}
			}
		}
		return false;
	}
	
	public boolean similar(Transcript t, double maxFraction)
	{
		//otherwise we check
		HashSet<Integer> borders = new HashSet<Integer>();
		int count=0;
		for(Exon e: exons)
		{
			if(count > 0)
				borders.add(e.start);
			if(count < exons.size())
				borders.add(e.end);
			count++;
		}
		count = 0;
		int coinciding = 0;
		int total = 0;
		for(Exon e: t.exons)
		{
			if(count > 0)
			{
				if(borders.contains(e.start))
					coinciding++;
				total++;
			}
			if(count < exons.size())
			{
				if(borders.contains(e.end))
					coinciding++;
				total++;
			}
			count++;
			borders.remove(e.start);
			borders.remove(e.end);
		}
		/*
		for(int b: borders)
		{
			if(b > t.getLeftMostEnd() && b < t.getRightMostStart())
				totalThis++; // the border was within the test transcript so this means the test transcript was missing an exon. 
		}
		*/
		double fractionTheSame = (1.0*(coinciding+2))/(total+2);
		if(fractionTheSame > maxFraction)
			return true;
		else
			return false; //enough different borders for this to quality as different
		
	}

	private int getRightMostStart() {
		// TODO Auto-generated method stub
		return Math.max(exons.getFirst().start,exons.getLast().start);
	}

	private int getLeftMostEnd() {
		// TODO Auto-generated method stub
		return Math.min(exons.getFirst().end,exons.getLast().end);
	}

	public double getSdIntensity() {

		double avgIntensity = getAvgIntensity();
		double sdIntensity = 0;
		for(Exon e: exons)
		{
			sdIntensity += Math.pow(e.getDensity()-avgIntensity,2);
		}
		return Math.sqrt(sdIntensity/exons.size())*TPM;
	}

	public int getExonCount() {
		return exons.size();
	}

	
	public String toString()
	{
		if(exons.size() == 0) return "";
		String result = "";

		for(Exon e: exons)
		{
			result+=e.start+"-"+e.end+" ";
		}
		return result;
	}

	public double getTPM() {
		return TPM;
	}
	
	public double getFPKM(){
		return FPKM;
	}

	public boolean contains(Exon current) {
		for(Exon e: exons)
		{
			if(e.start == current.start && e.end == current.end)
				return true;
		}
		return false;
	}

	public int length() {
		int result = 0;
		for(Exon e: exons)
			result+=(e.end-e.start);
		return result;
	}

	/**
	 * Returns the genomic position, taking splicing into account by traversing the exons.
	 * 
	 * @param offset: DNA sequence offset along transcript sequence
	 * @param startPad: added sequence pad from upstream transcript expansion
	 * @return
	 */
	public int getRealPos(int offset, int startPad, boolean plus) {
		
		if(offset <= startPad)
		{
			if(plus)
				return getStart()-startPad+offset;
			else return getEnd()+startPad-offset;
		}
		//offset is bigger than startPad
		offset-=startPad;

		int exonId = 0;
		while(exonId < exons.size() &&offset > 0)
		{
			offset -= (exons.get(exonId).end-exons.get(exonId).start);  //jump to end of exon
			exonId++;
		}
		exonId--;
		if(plus)
		{
			return exons.get(exonId).end+offset;
		}
		else
			return exons.get(exonId).start-offset;
	}

	public LinkedList<Exon> getSubset(int left, int right, int uPad, int dPad,boolean plus) {
		LinkedList<Exon> result = new LinkedList<Exon>();
		int index = 0;
		for(Exon e: exons)
		{
			if(e.overlaps(left, right))
			{
				//lets include the overlapping part
				int a = Math.max(left, e.start);
				int b = Math.min(right, e.end);
				if(a == b)
					continue;
				int size = b-a+1;
				if(index == 0 )  //first exon in transcript (left most in plus, rightmost in minus)
				{
					if(plus&& a == e.start)
						a-= uPad;
					else if(!plus&& b == e.end)
						b+=uPad;
				}
				if(index == exons.size()-1) //last exon
				{
					if(plus&& b == e.end)
						b+= dPad;
					else if(!plus&& a == e.start)
						a-=dPad;
				}
				Exon nexon = new Exon(a,b,e.plus,(e.count*size)/(e.end-e.start));
				result.add(nexon);
			}
			index++;
		}
		return result;
	}

	public double getExonLength() {
		double size = 0;
		for(Exon e: exons)
			size+=e.end-e.start;
		return size;
	}

	//reorders the exons
	public void flipExons() {
		LinkedList<Exon> temp = new LinkedList<Exon>();
		while(exons.size() > 0)
			temp.add(exons.removeLast());
		exons = temp;
	}

	public void trimExtendedReads(int ext) {
		Exon leftMost = null;
		Exon rightMost = null;
		for(Exon e: exons)
		{
			if(leftMost== null || leftMost.start > e.start)
				leftMost = e;
			if(rightMost == null || rightMost.end < e.end)
				rightMost = e;
		}
		leftMost.start=Math.min(leftMost.start+ext,leftMost.end-1);
		rightMost.end=Math.max(rightMost.end-ext,rightMost.start+1);
	}

	public int getMinExonLength() {
		int min = Integer.MAX_VALUE;
		for(Exon e: exons)
		{
			if(Math.abs(e.start-e.end) < min) min = Math.abs(e.start-e.end);
		}
		return min;
	}

	public double getMinExonIntensity() {
		double min = Double.POSITIVE_INFINITY;
		for(Exon e: exons)
		{
			if(e.getDensity() < min) min = e.getDensity();
		}
		return  min;
	}
}
