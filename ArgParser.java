/* Copyright (c)  2016   The Salk Institute for Biological Studies.
	All Rights Reserved
	 
	Permission to copy, modify and distribute any part of this MAPS for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
	Those desiring to incorporate this MAPS into commercial products or use for commercial purposes should contact the Technology Transfer Office, The Salk Institute for Biological Studies, La Jolla, 10010 N Torrey Pines Rd., La Jolla, CA 92037, Ph: (858) 453-4100.
	IN NO EVENT SHALL THE SALK INSTITUTE FOR BIOLOGICAL STUDIES BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS MAPS, EVEN IF THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	THE MAPS PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE SALK INSTITUTE FOR BIOLOGICAL STUDIES HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE SALK INSTITUTE FOR BIOLOGICAL STUDIES MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, 
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE MAPS WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/
import java.util.HashMap;
import java.util.LinkedList;

/**
 * Class designed for parsing arguments in a consistent manner. Register args, process a list of args, then ask for each arg in turn.
 * @author Maxim Shokhirev (C) 2016
 *
 */
public class ArgParser {
	private HashMap<String,String> args = new HashMap<String,String>();
	private LinkedList<String> list = new LinkedList<String>();
	private String description = "";
	private LinkedList<String> descriptions = new LinkedList<String>();
	
	public ArgParser(String description)
	{
		this.description = description;
	}
	
	public void registerArg(String name, String defaultValue, String description){
		args.put(name, defaultValue);
		descriptions.add("-"+name+"\t"+description+" (Default="+defaultValue+")");
	}
	
	
	public void printUsage()
	{
		System.err.println(description);
		System.err.println("\nOptions:");
		for(String str: descriptions)
		{
			System.err.println(str);
		}
		System.exit(1);
	}
	
	public void parseArgs(String[] args)
	{
		for(int i = 0; i < args.length; i++)
		{
			try{
			if(args[i].startsWith("-"))
			{
				String name = args[i].substring(1);

				if(i < args.length-1 && !args[i+1].startsWith("-"))
				{
					String value = args[i+1];
					i++;
					this.args.put(name,value);
				}
				else
					this.args.put(name,"true");
			}
			else
				list.add(args[i]);
			}catch(Exception e)
			{
				System.err.println("Error parsing: "+args[i]);
				printUsage();
			}
		}	
	}
	
	public String get(String key)
	{
		if(args.get(key)!= null)
			return args.get(key);
		return null;
	}
	
	public int getAsInt(String key)
	{
		if(args.get(key)!= null)
			return Integer.parseInt(args.get(key));
		return 0;
	}
	
	public double getAsDouble(String key)
	{
		if(args.get(key)!= null)
			return Double.parseDouble(args.get(key));
		return 0.0;
	}
	
	public boolean getAsBoolean(String key)
	{
		if(args.get(key)!= null)
			return Boolean.parseBoolean(args.get(key));
		return false;
	}
	
	public LinkedList<String> getList()
	{
		return list;
	}
	
	/**
	 * Prints values for all args: arg->value
	 * @return
	 */
	public LinkedList<String> getAllArgs()
	{
		LinkedList<String> result = new LinkedList<String>();
		for(String key: args.keySet())
		{
			result.add(key+"->"+args.get(key));
		}
		return result;
	}

}
