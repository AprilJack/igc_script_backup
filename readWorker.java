/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;


/**
 * Read workers are the workhorse of MAPS. They assemble reads into exons, exons into transcripts, 
 * filter exons/junctions/transcripts, split up exons by junctions, and produce a list of assembled filtered transcripts.
 * Stringency and overlap control the assembly and filtering. Designed on work on a specific set of reads coming from one 
 * strand of one chr such that they can be run in parallel to process an entire transcriptome simultaneously. 
 * @author Maxim Shokhirev. (C) 2015-2018
 *
 */
public class readWorker extends Thread {


	
	HashMap<String,LinkedList<int[]>> reads = null;
	
	HashMap<String,LinkedList<Exon>> exons  = null;
	HashMap<String,HashMap<Integer,Integer>> allJunctions = null;
	HashMap<String,HashMap<Integer,Integer>> junctions = new HashMap<String,HashMap<Integer,Integer>>();
	HashSet<Integer> spoints = new HashSet<Integer>();
	LinkedList<Exon> lastExons = new LinkedList<Exon>();
	LinkedList<Transcript> assembledTranscripts = new LinkedList<Transcript>();  //these are all of the transcripts that need to be printed out!
	int count = 0;
	boolean filter = false;
	boolean verbose = false;
	boolean stranded = true;
	String chrStr;
	int jfocus = 1;
	double stringency = 1;
	double mult = 1;
	double overlapMax = 0;
	int ext = 0;
	boolean keepReads = false;
	
	public readWorker(MAPS maps, String chrStr)
	{
		this.reads = maps.reads;
		this.mult = maps.mult;
		this.exons = maps.exons;
		this.allJunctions = maps.junctions;
		this.chrStr = chrStr;
		this.filter = maps.filter;
		this.verbose = maps.verbose;
		this.stringency = Math.min(0.99,Math.max(0.01,maps.stringency));
		this.jfocus=(int) (2+stringency*48);
		this.overlapMax = maps.overlapMax;
		this.ext = maps.ext;
		this.keepReads= maps.ap.getAsDouble("diversity")==1.0;
		maps.transcripts.put(chrStr,assembledTranscripts);
		this.start();
	}
	
	public void run()
	{
		try{
			
			boolean plus = true;
			if(chrStr.endsWith("-"))
				plus = false;
	
			HashMap<String,LinkedList<Exon>> map = exons;
			LinkedList<Exon> knownExons = new LinkedList<Exon>();
			map.put(chrStr,knownExons);
			LinkedList<int[]> theReads = reads.get(chrStr);
			LinkedList<Exon> lastExons = new LinkedList<Exon>();
			LinkedList<int[]> keptReads = new LinkedList<int[]>();
			while(theReads.size() > 0)
			{
				count++;
				int[] read = theReads.removeFirst();
				if(keepReads) keptReads.add(read);
				for(int r= 0; r < read.length; r+=2)
				{
					//lets find the exons that this read intersects with... since this list is sorted we know that we can just search the last few
					boolean found = false;
					LinkedList<Exon> tempLast = new LinkedList<Exon>();
					LinkedList<Exon> overlapping = new LinkedList<Exon>(); //the list of overlapping exons
					for(Exon e: lastExons)
					{
						if(e.overlaps(read[r], read[r+1],ext))  //check if one of
						{
							found = true;
							overlapping.add(e);
						}
						//reads are sorted by their first coord 
						if(r!=0 || e.end >= read[r]) //if the last exon is still relevant otherwise remove from list. Only remove if the first region starts afterward!
							tempLast.add(e);  
					}
					lastExons=tempLast;
					if(!found)
					{
						//adding a brand new exon
						Exon e = new Exon(read[r],read[r+1],plus,1);
						synchronized(map)
						{
							//add in order
							if(map.get(chrStr) == null)
							{
								map.put(chrStr,new LinkedList<Exon>());
								knownExons = map.get(chrStr);
							}
							knownExons.add(e);						
						}
						lastExons.add(e);
					}
					else
					{
						//we overlap with some previous exon(s)
						if(overlapping.size() == 1) //99% of the time
						{
							Exon overlap = overlapping.getFirst();
							overlap.end = Math.max(overlap.end, read[r+1]);
							overlap.start = Math.min(overlap.start, read[r]);
							overlap.count++;
						}
						else
						{
							//overlapping with multiple exons... lets tie them together (possible under rare circumstances)
							
							int min = Integer.MAX_VALUE;
							int max = Integer.MIN_VALUE;
							float total = 0;
							for(Exon overlap: overlapping)
							{
								if(overlap.start < min) min = overlap.start;
								if(overlap.end > max) max = overlap.end;
								total+=overlap.count;
								lastExons.remove(overlap);
								knownExons.remove(overlap);
							}
							min = Math.min(min,read[r]);
							max = Math.max(max,read[r+1]);
							Exon e = new Exon(min,max,plus,total+1);
							synchronized(map)
							{
								//add in order
								if(map.get(chrStr) == null)
								{
									map.put(chrStr,new LinkedList<Exon>());
									knownExons = map.get(chrStr);
								}
								knownExons.add(e);						
							}
							lastExons.add(e);
						}
					}
					if(r > 0)  //we are looking at a junctioned region
					{
						//the expansion of the exon and making new exons has already been handled in the part above
						synchronized(junctions)
						{
							if(read[r]-read[r-1] < 300000)  //the gap is small enough to be acceptable
							{  
								if(plus)
								{
									if(junctions.get(chrStr+read[r-1]) == null) junctions.put(chrStr+read[r-1],new HashMap<Integer,Integer>());
									if(junctions.get(chrStr+read[r-1]).get(read[r]) == null)
									{
										junctions.get(chrStr+read[r-1]).put(read[r],1);
						
									}
									else junctions.get(chrStr+read[r-1]).put(read[r],junctions.get(chrStr+read[r-1]).get(read[r])+1);
								}
								else
								{
									if(junctions.get(chrStr+read[r]) == null) junctions.put(chrStr+read[r],new HashMap<Integer,Integer>());
									if(junctions.get(chrStr+read[r]).get(read[r-1]) == null)
									{
										junctions.get(chrStr+read[r]).put(read[r-1],1);
								
									}
									else junctions.get(chrStr+read[r]).put(read[r-1],junctions.get(chrStr+read[r]).get(read[r-1])+1);
								}
								
							}
						}
					}
				}
			}//while have reads
			int filteredJunctions = 0;
			int totalJunctions = 0;
			int junctionPoints = 0;
			if(filter)
			{
				LinkedList<Exon> filteredExons = new LinkedList<Exon>();
				for(Exon e: knownExons)
				{
					if(!e.filter())
					{
						if(e.start <= e.end)
						{
							filteredExons.add(e);					
						}
						else
						{
							System.err.println("Exon ends before it starts: "+chrStr+" "+e.start+"-"+e.end);
						}
						
					}
				}
				map.put(chrStr, filteredExons);
				if(verbose)System.err.println("Filtered "+(knownExons.size()-filteredExons.size())+" poor quality exons from "+chrStr+" total exons now: "+filteredExons.size());
				//now lets filter the junctions
				HashMap<String,HashMap<Integer,Integer>> junctions2 = new HashMap<String,HashMap<Integer,Integer>>();
				for(String str: junctions.keySet())
				{
					totalJunctions+=junctions.get(str).size();  //find the total number of end points for this junction
					String fromStr=str.split("[+-]")[1];
					int from = Integer.parseInt(fromStr);   
					String chr = str.substring(0,str.length()-fromStr.length());
					double ourBest = 0;
					for(int where: junctions.get(str).keySet())
					{
						double count = (1.0*(junctions.get(str).get(where)))/Math.sqrt(Math.abs(from-where));  //finding the most supported end-point
						if(count > ourBest)
						{
							ourBest =count;
						}
 
					}
					//ourBest is the highest count for this junction
					boolean neighbor_is_better = false;
					for(int i =-jfocus; i <=jfocus; i++)   //lets check if there is a better junc starting near our start
					{
						if(i != 0 &&junctions.get(chr+(from+i))!= null)
						{
							for(int where: junctions.get(chr+(from+i)).keySet())
							{
								if(junctions.get(chr+(from+i)).get(where)/Math.sqrt(Math.abs(from-where)) > ourBest)  //check if the neighbor has a better end point
								{
									neighbor_is_better =true;
									break;
								}
							}
						}
						if(neighbor_is_better)
							break;
					}
					if(!neighbor_is_better && ourBest > 0.005*mult*stringency*stringency)  //if the junction is the best in the vicinity AND it has at least a few counts
					{
						HashMap<Integer,Integer> temp2 = new HashMap<Integer,Integer>();
						double bestScore = 0;
						for(int where: junctions.get(str).keySet()) //for each end point find the best count
						{				
							double current = junctions.get(str).get(where)/Math.sqrt(Math.abs(from-where));
							if(current > bestScore)
							{
								bestScore = current;
							}
						}
						for(int where: junctions.get(str).keySet())
						{
							double current = junctions.get(str).get(where)/Math.sqrt(Math.abs(from-where));
							if(current >= stringency*bestScore|| current == bestScore)   //save the junction with all "decent" end points
							{
								temp2.put(where,junctions.get(str).get(where));
							}
						}
						if(temp2.size() > 0)junctions2.put(str, temp2);
						filteredJunctions+=junctions.get(str).size()-temp2.size();
					}
					else
						filteredJunctions+=junctions.get(str).size();
				}
				junctions=junctions2;
				//filterJunctionsDownTo(2);
				
			}
			allJunctions.putAll(junctions);
			if(verbose)System.err.println(String.format("%s had %d junctions of which %3.1f%% were filtered.",chrStr,totalJunctions,((100.0*filteredJunctions)/totalJunctions)));
			Collections.sort(exons.get(chrStr)); // sort the exons based on their left-most or right-most coord (makes sense)
			//now lets break up the exons based on splice junctions
			LinkedList<Exon> elist = exons.get(chrStr);
			int startCount = elist.size();
			LinkedList<Exon> elist2 = new LinkedList<Exon>();
			if(verbose)System.err.println("Splitting up "+elist.size()+" exons on "+chrStr);
			HashSet<Integer> endPoints = new HashSet<Integer>();
			while(elist.size() > 0)
			{
				Exon e = elist.removeFirst();
				boolean added = false;
				if(plus)
				{
					for(int splice = e.start; splice <= e.end; splice++)
					{
						if(junctions.get(chrStr+splice) != null)
						{
							//lets see where it goes
							for(int where: junctions.get(chrStr+splice).keySet())
							{

								if(where > e.end)  //goes to outside the current exon so split this one up
								{
									//first and last part added to the list
									int size1 =splice-e.start;
									int size2 =e.end-splice;
									if(size1 > 2)
										elist2.add(new Exon(e.start,splice,e.plus,e.count*(size1/((float)e.end-e.start))));
									if(size2 > 2)
										elist2.add(new Exon(splice,e.end,e.plus,e.count*(size2/((float)e.end-e.start))));
									//lets add the end point for future cutting
									endPoints.add(where);
									added=true;
									junctionPoints+=2;
		
								}
							}
						}
					}
				}
				else // minus
				{
					for(int splice = e.end; splice >= e.start; splice--)
					{
						if(junctions.get(chrStr+splice) != null)
						{
							//lets see where it goes
							for(int where: junctions.get(chrStr+splice).keySet())
							{
								if(where < e.start)  //goes to outside the current exon so split this one up
								{
									//first and last part added to the list
									int size1 =e.end-splice;   //after the splice
									int size2 =splice-e.start;  //before the splice
									if(size1 > 2)
										elist2.add(new Exon(splice,e.end,e.plus,e.count*(size1/((float)e.end-e.start))));
									if(size2 > 2)
										elist2.add(new Exon(e.start,splice,e.plus,e.count*(size2/((float)e.end-e.start))));
									//lets add the end point for future cutting
									endPoints.add(where);
									added=true;
									junctionPoints+=2;
								}
							}
						}
					}
				}
				if(!added)
					elist2.add(e); //in case we have internal splice junctions and end points that need to be handled
			}
			//Collections.sort(elist2); don't need to sort because the chop order doesn't matter. only the final products
			HashSet<String> added = new HashSet<String>();
			while(elist2.size() > 0)
			{
				Exon e = elist2.removeFirst();
				boolean spliced = false;
				if(plus)
				{
					for(int splice = e.start; splice < e.end; splice++)
					{
						if(junctions.get(chrStr+splice) != null)
						{
							//lets see where it goes
							for(int where: junctions.get(chrStr+splice).keySet())
							{
								if(where > e.start && where < e.end) //going to inside this exon
								{
									//first part added to new list
									Exon addBack = new Exon(where,e.end,e.plus,e.count*((e.end-where)/((float)e.end-e.start)));
									elist.add(new Exon(e.start,splice,e.plus,e.count*((splice-e.start)/((float)e.end-e.start))));
									elist.add(new Exon(splice,where,e.plus,e.count*((where-splice)/((float)e.end-e.start))));
									elist.add(e);
									if(!added.contains(addBack.start+"-"+addBack.end))
									{
										elist2.addFirst(addBack);
										added.add(addBack.start+"-"+addBack.end);
									}
									spliced = true;
									junctionPoints+=2;
								}
							}
						}
	 	 			}
				}
				else
				{
					for(int splice = e.end; splice > e.start; splice--)
					{
						if(junctions.get(chrStr+splice) != null)
						{
							//lets see where it goes
							for(int where: junctions.get(chrStr+splice).keySet())
							{
								if(where > e.start && where < e.end)  //going inside
								{
									//first part added to new list
									Exon addBack = new Exon(e.start,where,e.plus,e.count*((where-e.start)/((float)e.end-e.start)));
									elist.add(new Exon(where,splice,e.plus,e.count*((splice-where)/((float)e.end-e.start))));
									elist.add(new Exon(splice,e.end,e.plus,e.count*((e.end-splice)/((float)e.end-e.start))));
									elist.add(e);  //add the full length again in case we want to skip this junction completely
									if(!added.contains(addBack.start+"-"+addBack.end))
									{
										elist2.addFirst(addBack);
										added.add(addBack.start+"-"+addBack.end);
									}
									spliced = true;
									junctionPoints+=2;
								}
							}
						}
	 	 			}
				}
				if(!spliced)
					elist.add(e); //for final endpoint cutting
				
			}
			added.clear();
			added = null;
			Collections.sort(elist);
			//at this point elist2 is empty
			//don't need to sort since we are not adding anything back in.
			while(elist.size() >0)
			{
				Exon e = elist.removeFirst();
				//boolean spliced = false;
				for(int splice = e.start; splice <= e.end; splice++)
				{
					if(endPoints.contains(splice))
					{
						
						int size1 = splice-e.start;
						int size2 = e.end-splice;
						if(size1 > 2)
							elist2.add(new Exon(e.start,splice,e.plus,e.count*(size1/((float)e.end-e.start))));
						if(size2 > 2)
							elist2.add(new Exon(splice,e.end,e.plus,e.count*(size2/((float)e.end-e.start))));
						//spliced = true;
					}
				}
				//if(!spliced)
					elist2.add(e);
			}
			//at this point elist is empty and elise2 contains the split up exons
			elist = new LinkedList<Exon>();
			HashSet<String> startsEnds = new HashSet<String>();
			for(Exon e: elist2)
			{
				if(!startsEnds.contains(e.start+";"+e.end))
				{
					startsEnds.add(e.start+";"+e.end);
					elist.add(e);
				}
			}
			elist2.clear();
			for(Exon e: elist)
			{
				if(e.start < e.end)
					elist2.add(e);
			}
			Collections.sort(elist2);
			exons.put(chrStr, elist2);
			if(verbose)System.err.println("Split up "+startCount+" exons on "+chrStr+" into "+elist2.size()+" along "+junctionPoints+" junction points.");
			//Now lets generate recursive transcripts for outputing
			LinkedList<Exon> tmpCopy = new LinkedList<Exon>();
			tmpCopy.addAll(elist2);
			LinkedList<Transcript> lastTranscripts = new LinkedList<Transcript>();
			int multi_exon=0 , stringent = 0, unique = 0, branches = 0,bcount=0;
			Exon lastExon = null;
			while(tmpCopy.size() > 0)
			{
				try{
					Exon current = tmpCopy.removeFirst();  //for each exon
					if(lastExon != null && ((plus && current.start == lastExon.start)||(!plus && current.end == lastExon.end)))
						continue;
					lastExon = current;
					LinkedList<Exon> allExons = new LinkedList<Exon>();
					allExons.addAll(tmpCopy);				
						
					LinkedList<Transcript >transcripts = new LinkedList<Transcript>();
					HashSet<String> tried = new HashSet<String>();
					while(true)
					{
						Transcript thusFar = new Transcript();
						long tic = System.currentTimeMillis();
						thusFar = getBestTranscript(plus, chrStr, current, allExons, thusFar, tried, "");
						long toc = System.currentTimeMillis()-tic;
						transcripts.add(thusFar);
						if(toc > 1000)
						{
							if(verbose)System.err.println("Encountered a very complicated transcript. Moving on in the interest of time...");
							break;
						}
 						if(thusFar.getExonCount() == 1 || transcripts.size() > 1000)
							break;
					}
					Collections.sort(transcripts,new Comparator<Transcript>(){

						@Override
						public int compare(Transcript o1, Transcript o2) {
							// TODO Auto-generated method stub
							double c1 = o1.getExonCount();
							double c2 = o2.getExonCount();
							if(c1 < c2) return 1;
							else if(c2 > c2) return -1;
							else return 0;
						}
						
					});
					//int truncatedTranscriptSize = transcripts.size();
					/*
					if(transcripts.size() > 1000)
					{
						LinkedList<Transcript> survivors = new LinkedList<Transcript>();
						while(survivors.size() < Math.min(truncatedTranscriptSize,Math.max(1000,Math.sqrt(truncatedTranscriptSize))))
							survivors.add(transcripts.removeFirst());
						transcripts = survivors;
					}
					*/
					branches += transcripts.size();
					bcount++;
					//lets sort the transcripts by exonCount
					//if(verbose) System.err.println(chrStr+" had a total of "+origTranscriptSize+" hypothetical transcripts of which the top "+transcripts.size()+" were kept for further filtering.");
					for(Transcript transcript: transcripts)
					{
						//end of exon search loop
						//lets output the transcript if it is of decent quality

						if(transcript.getMinExonLength() > 4 && (transcript.getExonCount()> 2 || transcript.getAvgIntensity()*Math.sqrt(transcript.getExonLength()+1)>0.005*mult*stringency))
								
						{
							if(transcript.getExonCount() > 1)
								multi_exon++;
							stringent++;
							//lets check for overlap to previous transcripts and only keep novel ones
							boolean good = true;
							LinkedList<Transcript> remaining = new LinkedList<Transcript>();

							
							for(Transcript old: lastTranscripts)
							{
								//if the old transcript covers the new transcript
								
								if((transcript.getExonCount()==1&& old.exonsCover(transcript))||
										(transcript.getExonCount() > 1 &&(old.similar(transcript,1.0-0.25*stringency)))) 
								{
									good = false;
									break;
								}
							}

							
							for(Transcript old: lastTranscripts)
							{
								
								if((plus && old.getEnd() > transcript.getStart())||(!plus && old.getStart() < transcript.getEnd()))
									remaining.add(old); //keep it around since it could still be overlapping for a new transcript
							}
														
							lastTranscripts = remaining;			
							//otherwise lets add this guy into the list
							if(good) 
							{
								lastTranscripts.addFirst(transcript);
								assembledTranscripts.add(transcript);
								unique++;
							}
						}
					}
					count+=transcripts.size();
				}
				catch(Exception e){
					System.err.println("Encountered an error during transcript assembly: "+e.getMessage()+" continuing bravely!");
				}
			}
			if(verbose)System.err.println(String.format("%s had %3.3f%% unique %3.3f%% stringent %3.3f%% multi-exon transcripts out of a total of %d (%d remaining). Avg transcripts per bundle: %3.3f",chrStr,(100.0*unique)/count,(100.0*stringent)/count,(100.0*multi_exon)/count,count,assembledTranscripts.size(),(1.0*branches)/bcount));
			reads.put(chrStr, keptReads);
			
		}catch(Exception e){
			e.printStackTrace();
		}
		//finding abundances
		computeTranscriptAbundances();

	}
	/*
	private void filterJunctionsDownTo(int maxBranches)
	{
		HashMap<String, HashMap<Integer,Integer>> filtered = new HashMap<String,HashMap<Integer,Integer>>();
		int removed = 0;
		for(String fromStr: junctions.keySet())
		{
			ArrayList<Integer> toList = new ArrayList<Integer>();
			ArrayList<Integer> toCount = new ArrayList<Integer>();
			for(int to: junctions.get(fromStr).keySet())
			{
				int i = 0;
				for(i=0; i < toCount.size(); i++)
				{
					if(toCount.get(i) < junctions.get(fromStr).get(to))
						break;
				}
				toList.add(i,to);
				toCount.add(i,junctions.get(fromStr).get(to));
			}
			filtered.put(fromStr,new HashMap<Integer,Integer>());
			for(int i = 0; i < Math.min(maxBranches,toList.size()); i++){
				filtered.get(fromStr).put(toList.get(i), toCount.get(i));
			}
			removed+=Math.max(0,toList.size()-maxBranches);
		}
		System.err.println("Removed "+removed);
		junctions = filtered;
	}
	
	private int estimateGoodBranchPoints(Exon current,
			LinkedList<Exon> allExons, boolean plus, String chr, int branchesThusFar) throws OutOfMemoryError {

		int result = 1;
		if(branchesThusFar > 100000)
			return 1;  //just stop it early to avoid a complete explosion
		//LETS TRY TO FIND A JUNCTION STARTING AT THE END OF THIS EXON
		if(plus)  
		{
			HashMap<Integer,Integer> temp = junctions.get(chr+current.end);
			if(temp != null)
			{
				int max = 0;
				int sum = 0;
				for(int whereTo: temp.keySet())  //possible destinations from the end of this exon
				{
					if(temp.get(whereTo)> max)
						max = temp.get(whereTo);
					sum+=temp.get(whereTo);
				}
				for(int whereTo : temp.keySet()) //for each possible jump
				{
					double nextProb = ((double)temp.get(whereTo))/sum;  //probability of the jump given all support for all junctions
					if(nextProb >= Math.max(0.1,squaredStringency) || temp.get(whereTo)==max)//cutoff for low abundance (e.g. don't follow bad junction ends)
					{
						//lets recurse to the exon at the end of this junction
						ArrayList<Exon> candidates = new ArrayList<Exon>();
						ArrayList<Double> scores = new ArrayList<Double>();
						double bestScore = 0;
						for(Exon e: allExons)
						{
							if(e.start == whereTo)//this exon starts where the junction leads!
							{
								double best = 0;
								if(junctions.get(chr+e.end) != null) //lets find the best exon starting at the end of the junction
								{ //best is defined as the exon with the highest support for a junction starting at its end

									for(int where: junctions.get(chr+e.end).keySet())
									{
										if((1.0*junctions.get(chr+e.end).get(where))/(e.end-e.start) > best)
											best=(1.0*junctions.get(chr+e.end).get(where))/(e.end-e.start);
									}

								}
								if(candidates.size() < 3)
								{
									candidates.add(e);  //add the exon to candidate list regardless using their best jump score as the score
									scores.add(best);
									if(best > bestScore)
										bestScore = best;
								}
							}
							else if(e.start > whereTo)
								break;  //all the rest will also be more
						}
						for(int i = 0; i < candidates.size(); i++)
						{
							if(scores.get(i)== bestScore)  //for all exons starting at the end of this junction, if they have support for their next jump or if this is the ending exon
							{
								result+=candidates.size()*estimateGoodBranchPoints(candidates.get(i),allExons,plus,chr,branchesThusFar*candidates.size());
								if(result == 0)
									break;
							}
						}
					}					
				}
			}
		}
		else
		{
			HashMap<Integer,Integer> temp = junctions.get(chr+current.start);
			if(temp != null)
			{
				int max = 0;
				int sum = 0;
				for(int whereTo: temp.keySet())
				{
					if(temp.get(whereTo)> max)
						max = temp.get(whereTo);
					sum+=temp.get(whereTo);
				}
				for(int whereTo : temp.keySet())
				{
					double nextProb = ((double)temp.get(whereTo))/sum;
					if(nextProb >= stringency || temp.get(whereTo)==max)//cutoff for low abundance (e.g. don't follow bad junction ends)				
					{
						//lets recurse to the exon at the end of this junction
						ArrayList<Exon> candidates = new ArrayList<Exon>();
						ArrayList<Double> scores = new ArrayList<Double>();
						double bestScore = 0;
						for(Exon e: allExons)
						{
							if(e.end == whereTo)//this exon starts where the junction leads!
							{
								double best = 0;
								if(junctions.get(chr+e.start) != null)
								{

									for(int where: junctions.get(chr+e.start).keySet())
									{
										if((1.0*junctions.get(chr+e.start).get(where))/(e.end-e.start) > best)
											best=(1.0*junctions.get(chr+e.start).get(where))/(e.end-e.start);
									}

								}
								if(candidates.size() < 2)
								{
									candidates.add(e);
									scores.add(best);
									if(best > bestScore)
										bestScore = best;
								}
							}
							else if(e.end < whereTo)//assume no exon larger than 100k in length exists
								break;  //all the rest will also be less (poor man's alternative to alternate structure for exons on negative strand)
						}
						for(int i = 0; i < candidates.size(); i++)
						{
							if(scores.get(i) == bestScore)  //for all exons starting at the end of this junction, if they have support for their next jump or if this is the ending exon
							{
								result+=candidates.size()*estimateGoodBranchPoints(candidates.get(i),allExons,plus,chr,branchesThusFar*candidates.size());
								if(result == 0)
									break;
							}
						}
					}
				}
			}
		}
		if(result > 100000000 || result == 0)
			return 0;
		return result;
	}
	*/
	/**
	 * Goes through each chromosome strand and finds clusters of overlapping transcripts. Those transcripts are then divided into
	 * virtual exons and the true transcript abudnance is calculated based on pre-recorded exon read counts. This is done
	 * using a simulated annealing optimization on all transcripts in the cluster. The simplification is that reads are on average evenly
	 * distributed across exons. This may not be true in general so take with a grain of salt.
	 */
	private void computeTranscriptAbundances() {
		double totalFractions = 0;
		//now lets go through and figure out boundaries for all exons in all of these transcripts
		HashSet<Integer> points = new HashSet<Integer>();
		for(Transcript t:assembledTranscripts)
		{
			for(Exon e: t.getExons())
			{
				points.add(e.start);
				points.add(e.end);
			}
		}
		HashMap<String,LinkedList<Transcript>> vexons = new HashMap<String,LinkedList<Transcript>>(); //START,END pair
		HashMap<String,Double> vexonCount = new HashMap<String,Double>();
		LinkedList<Transcript> ts = new LinkedList<Transcript>();
		//ok now lets go through each transcript and populate a set of virtual exons until there is no overlap.
		//then we want to do the optimization on the current cluster of virtual-exons and transcripts
		int count = 0;
		for(Transcript t: assembledTranscripts)
		{
			HashMap<String,LinkedList<Transcript>> vexons2 = new HashMap<String,LinkedList<Transcript>>(); //START,END pair
			HashMap<String,Double> vexonCount2 = new HashMap<String,Double>();
			count++;
			boolean overlaped = false;
			for(Exon e: t.getExons())
			{
				int left = -1;
				//lets divide it up into virtual exons
				for(int x = e.start; x <= e.end; x++)
				{
					if(points.contains(x))
					{
						if(left != -1)
						{
							String vexon = left+","+x;
							vexons2.put(vexon,new LinkedList<Transcript>());
							vexons2.get(vexon).add(t);
							vexonCount2.put(vexon,(double) (e.count)/(e.end-e.start));
							if(vexons.get(vexon) != null)
								overlaped =true;
						}
						left = x;
					}
				}
			}
			//now lets check if we had any overlap with previous transcripts in the cluster
			if(overlaped)
			{
				for(String vexon: vexons2.keySet())
				{
					if(vexons.get(vexon) == null) {
						vexons.put(vexon,new LinkedList<Transcript>());
						vexonCount.put(vexon,vexonCount2.get(vexon));
					}
					vexons.get(vexon).add(t);
				}

			}
			//else lets do the optimization and replace swap out the maps
			else
			{
				if(vexons.size() > 0) //first time around we wont have anything
				{
					//lets run simulated annealing on these clusters of transcripts
					double score =runSimulatedAnnealing(ts,vexons,vexonCount);
					//System.err.println((100.0*count)/assembledTranscripts.size()+" "+score);
					for(Transcript tr: ts)
						totalFractions+=tr.getFPKM();
					//if(verbose)System.err.println("Simulated annealing score was: "+score);
				}
				vexons = vexons2;
				vexonCount = vexonCount2;
				ts.clear();
			}
			ts.add(t);
		}
		double avgError = 0;
		int tcount = 0;

		for(Transcript t: assembledTranscripts)
		{
			t.setTPM((t.getFPKM()*1000000)/totalFractions);
			avgError+=t.getError();
			tcount++;
		}
		if(verbose)System.err.println(String.format("Finished assembling and quantifying %d transcripts on %s. Average quality: %3.3f",tcount,chrStr,avgError/tcount));
	}

	/**
	 * Will run simulated annealing optimization to find the optimal transcript count
	 * that results in the given virutal exon counts
	 * @param ts
	 * @param vexons
	 * @param vexonCount
	 * @return normalized squared error
	 */
	private double runSimulatedAnnealing(LinkedList<Transcript> ts,
			HashMap<String, LinkedList<Transcript>> vexons,
			HashMap<String, Double> vexonCount) {
		double T = 10.0;
		double score = 99999999;
		double totalCount = 0;
		for(double count: vexonCount.values())
			totalCount+=count;
		//initialize counts
		for(Transcript t: ts)
		{
			t.setFPKM(Math.random()*5);
		}
		int count = 10000;
		if(ts.size() ==1) count = 5000;
		double Tfactor = Math.exp(Math.log(0.001)/(Math.min(1000000,count*ts.size()))); //increasing to 100k improved the score by just 2%.
		long tic = System.nanoTime();
		while(T > 0.01)
		{
			//lets generate a neighbor 
			HashMap<Transcript, Double> oldScores = new HashMap<Transcript,Double>();
			for(Transcript t: ts)
			{
				oldScores.put(t,t.getFPKM());
				t.setFPKM(Math.max(0,t.getFPKM()*Math.random()*2+(Math.random()-Math.random())*T));
			}
			double currentScore = 0;
			for(String vexon: vexons.keySet())
			{
				double sum = 0;
				for(Transcript t: vexons.get(vexon))
				{
					sum+=t.getFPKM();
				}
				currentScore+=Math.pow((sum-vexonCount.get(vexon))/totalCount,2);
			}
			currentScore/=vexons.size();
			if(currentScore < score || Math.random() < Math.exp((score-currentScore)/(T*T)))
			{
				//keep the values
				score = currentScore;
			}
			else
			{
				//revert back
				for(Transcript t: oldScores.keySet())
				{
					t.setFPKM(oldScores.get(t));
					t.setError(-Math.log(score));
				}
			}
//			System.err.println(count+"\t"+T+"\t"+score);
			count++;
			T=T*Tfactor;
			if(System.nanoTime()-tic >10000000)
				Tfactor/=2;
			tic = System.nanoTime();
		}

		//System.err.println("T count="+ts.size()+" score was "+score);
		return score;
		
	}
	
	/*
	private class Cluster
	{
		LinkedList<Exon> exons = new LinkedList<Exon>();
		HashMap<Integer,HashMap<Integer,Integer>> junctions = new HashMap<Integer,HashMap<Integer,Integer>>();
		
		public Cluster()
		{
			//do nothing
		}
		
		public void addExon(Exon e)
		{
			exons.add(e);
		}
		
		public void addJunction(int from, int to, int reads)
		{
			if(junctions.get(from) == null)
				junctions.put(from,new HashMap<Integer,Integer>());
			junctions.get(from).put(to, reads);
		}
		
		
		 // Removes all but the top 2 junctions for any one starting position

		private void filterJunctions()
		{
			HashMap<Integer, HashMap<Integer,Integer>> filtered = new HashMap<Integer,HashMap<Integer,Integer>>();
			for(Integer start: filtered.keySet())
			{
				ArrayList<Integer> toList = new ArrayList<Integer>();
				ArrayList<Integer> toCount = new ArrayList<Integer>();
				for(int to: filtered.get(start).keySet())
				{
					int i = 0;
					for(i=0; i < toCount.size(); i++)
						if(toCount.get(i) < filtered.get(start).get(to))
							break;
					toList.add(i,to);
					toCount.add(i,filtered.get(start).get(to));
				}
				filtered.put(start,new HashMap<Integer,Integer>());
				for(int i = 0; i < Math.min(2,toList.size()); i++){
					filtered.get(start).put(toList.get(i), toCount.get(i));
				}
			}
			junctions = filtered;
		}
		
		
		 // Returns a list of transcripts

		public LinkedList<Transcript> getTranscripts(String chr, boolean plus)
		{
			filterJunctions();
			Collections.sort(exons);
			LinkedList<Transcript> transcripts = new LinkedList<Transcript>();
			for(Exon e: exons)
			{
				transcripts.addAll(recursiveTranscript(new Transcript(), e, exons, plus, chr, stringency));				
			}
			return transcripts;
		}
	}
	*/
	
	private Transcript getBestTranscript(boolean plus, String chr, Exon current, LinkedList<Exon> allExons, Transcript thusFar, HashSet<String> tried, String choices) {
		//LETS TRY TO FIND A JUNCTION STARTING AT THE END OF THIS EXON
		Transcript updatedThusFar = thusFar.duplicate();
		updatedThusFar.add(current);
		tried.add(choices);
		HashMap<Integer,Integer> temp = junctions.get(chr+current.end);
		if(!plus) temp = junctions.get(chr+current.start);
		if(temp != null)
		{
			int max = 0;
			int sum = 0;
			ArrayList<Integer> toList = new ArrayList<Integer>();
			ArrayList<Double> toScore = new ArrayList<Double>();
			for(int whereTo: temp.keySet())  //possible destinations from the end of this exon
			{
				if(temp.get(whereTo)> max)
					max = temp.get(whereTo);
				sum+=temp.get(whereTo);
			}
			for(int whereTo: temp.keySet())  //possible destinations from the end of this exon
			{
				double prob = ((double)temp.get(whereTo))/(sum*Math.sqrt(Math.min(Math.abs(whereTo-current.end),Math.abs(whereTo-current.start))));
				int index = 0; 
				while(index < toList.size())
				{
					if(toScore.get(index) < prob)
						break;
					else index++;
				}
				toList.add(index,whereTo);
				toScore.add(index,prob);
			}
			
			for(int i = 0; i < Math.min(toList.size(),3); i++) //for each possible jump
			{
				int whereTo = toList.get(i);
				double nextProb = toScore.get(i);  //probability of the jump given all support for all junctions
				if(nextProb >= stringency*stringency || temp.get(whereTo)==max)//cutoff for low abundance (e.g. don't follow bad junction ends)
				{
					//lets recurse to the exon at the end of this junction
					ArrayList<Exon> candidates = new ArrayList<Exon>();
					ArrayList<Double> scores = new ArrayList<Double>();
					double bestScore = 0;
					for(Exon e: allExons)
					{
						if((plus && e.start == whereTo)||(!plus && e.end == whereTo))//this exon starts where the junction leads!
						{
							double score = e.getDensity()/Math.abs(e.end-e.start);
							int index = 0;
							while(index < candidates.size())
							{
								if(scores.get(index) < score)
									break;
								else index++;
							}
							candidates.add(index,e);  //add the exon to candidate list regardless using their best jump score as the score
							scores.add(index,score);
							if(score > bestScore)
								bestScore = score;	
						}
						else if((plus &&e.start > whereTo)||(!plus && e.end < whereTo))
							break;  //all the rest will also be more
					}
					if(updatedThusFar.getExonCount() < 500)  //Titin has 364 exons
					{
						//lets check if we picked this combo before 
						
						for(int j = 0; j < Math.min(3,candidates.size()); j++)
						{
							if(!tried.contains(choices+"|"+j))
							{
								Transcript best = getBestTranscript(plus,chr,candidates.get(j),allExons,updatedThusFar,tried,choices+"|"+j);
								return best;
							}
						}
						return updatedThusFar; //we have hit a dead end in this tree so this is our finish case
					}
				}					
			}
		}
	

		return updatedThusFar;
	}
	
	
	/**
	 * Recursively assembles transcripts keeping stringency in mind. Combinatorial explosion is avoided by temporarily increasing stringency.
	 * @param current the current exon being handled
	 * @param allExons all exons on chr
	 * @param currentExons list of exons in order for this transcript
	 * @param intensities list of exon intensities in order for this transcript
	 * @param transcriptSeq  the transcript sequence in ACTG
	 * @param transcriptStart the start of this transcript in chr coords
	 * @param transcriptEnd the end of this transcript in chr coords
	 * @param plus directionality of transcript
	 * @param chr 
	 */
	/*
	private LinkedList<Transcript> recursiveTranscript(Transcript thusFar, Exon current, LinkedList<Exon> allExons,
			boolean plus,String chr, double currentStringency) throws OutOfMemoryError {
		thusFar.add(current);
		LinkedList<Transcript> result = new LinkedList<Transcript>();

		//LETS TRY TO FIND A JUNCTION STARTING AT THE END OF THIS EXON
		if(plus)  
		{
			HashMap<Integer,Integer> temp = junctions.get(chr+current.end);
			if(temp != null)
			{
				int max = 0;
				int sum = 0;
				for(int whereTo: temp.keySet())  //possible destinations from the end of this exon
				{
					if(temp.get(whereTo)> max)
						max = temp.get(whereTo);
					sum+=temp.get(whereTo);
				}
				for(int whereTo : temp.keySet()) //for each possible jump
				{
					double nextProb = ((double)temp.get(whereTo))/sum;  //probability of the jump given all support for all junctions
					if(nextProb >= currentStringency || temp.get(whereTo)==max)//cutoff for low abundance (e.g. don't follow bad junction ends)
					{
						//lets recurse to the exon at the end of this junction
						ArrayList<Exon> candidates = new ArrayList<Exon>();
						ArrayList<Double> scores = new ArrayList<Double>();
						double bestScore = 0;
						for(Exon e: allExons)
						{
							if(e.start == whereTo)//this exon starts where the junction leads!
							{
								double best = 0;
								if(junctions.get(chr+e.end) != null) //lets find the best exon starting at the end of the junction
								{ //best is defined as the exon with the highest support for a junction starting at its end

									for(int where: junctions.get(chr+e.end).keySet())
									{
										if((1.0*junctions.get(chr+e.end).get(where))/(e.end-e.start) > best)
											best=(1.0*junctions.get(chr+e.end).get(where))/(e.end-e.start);
									}

								}
								if(candidates.size() < 3)
								{
									candidates.add(e);  //add the exon to candidate list regardless using their best jump score as the score
									scores.add(best);
									if(best > bestScore)
										bestScore = best;	
								}
							}
							else if(e.start > whereTo)
								break;  //all the rest will also be more
						}
						if(thusFar.getExonCount() < 500)  //Titin has 364 exons
						{
							for(int i = 0; i < Math.min(3,candidates.size()); i++)
							{
								if(scores.get(i) >= bestScore*currentStringency || bestScore==0)  //for all exons starting at the end of this junction, if they have support for their next jump or if this is the ending exon
								{
									result.addAll(recursiveTranscript(thusFar.duplicate(),candidates.get(i),allExons,plus,chr,currentStringency));
								}
							}
						}
					}					
				}
			}
		}
		else
		{
			HashMap<Integer,Integer> temp = junctions.get(chr+current.start);
			if(temp != null)
			{
				int max = 0;
				int sum = 0;
				for(int whereTo: temp.keySet())
				{
					if(temp.get(whereTo)> max)
						max = temp.get(whereTo);
					sum+=temp.get(whereTo);
				}
				for(int whereTo : temp.keySet())
				{
					double nextProb = ((double)temp.get(whereTo))/sum;
					if(nextProb >= currentStringency || temp.get(whereTo)==max)//cutoff for low abundance
					{
						//lets recurse to the exon at the end of this junction
						ArrayList<Exon> candidates = new ArrayList<Exon>();
						ArrayList<Double> scores = new ArrayList<Double>();
						double bestScore = 0;
						for(Exon e: allExons)
						{
							if(e.end == whereTo)//this exon starts where the junction leads!
							{
								double best = 0;
								if(junctions.get(chr+e.start) != null)
								{

									for(int where: junctions.get(chr+e.start).keySet())
									{
										if((1.0*junctions.get(chr+e.start).get(where))/(e.end-e.start) > best)
											best=(1.0*junctions.get(chr+e.start).get(where))/(e.end-e.start);
									}

								}
								if(candidates.size() < 3)
								{
									candidates.add(e);
									scores.add(best);
									if(best > bestScore)
										bestScore = best;
								}

							}
							else if(e.end < whereTo)//negative strand exons are sorted by their right-most coord
								break;  //all the rest will also be less (poor man's alternative to alternate structure for exons on negative strand)
						}
						if(thusFar.getExonCount() < 500)
						{
							for(int i = 0; i < Math.min(3,candidates.size()); i++)
							{
								if(scores.get(i)>= bestScore*currentStringency || bestScore == 0)
								{
									result.addAll(recursiveTranscript(thusFar.duplicate(),candidates.get(i),allExons,plus,chr,currentStringency));
								}
							}
						}
					}
				}
			}
		}
		if(result.size() == 0)
			result.add(thusFar);
		return result;
	}
	*/
}
