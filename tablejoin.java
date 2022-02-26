

import java.awt.Component;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.ListCellRenderer;
import javax.swing.ListSelectionModel;

public class tablejoin 
{
	
	public tablejoin(String[] args)
	{
		if(args.length < 2)
		{
			System.err.println("Usage: tablejoin File1.csv File2.tab [...]");
			System.err.println("Combines rows starting with the same name in column 1 of each file.");
			System.err.println("Only rows with ids in all files are printed to stdout.");
			System.exit(1);
		}
		//load in the first table right away
		LinkedList<String[]> table = loadTable(args[0]);
		System.err.println("First table had "+table.size());
		HashMap<String,Integer> idmap = new HashMap<String,Integer>();
		for(int i = 0; i < table.size(); i++)
			idmap.put(table.get(i)[0], i);
		int cols = table.get(1).length;
		for(int i = 1; i < args.length; i++)
		{
			LinkedList<String[]> table2 = loadTable(args[1]);
			System.err.println("Table "+(i+1)+" had "+table2.size());
			cols += table2.get(1).length-1;
			for(int j = 0; j < table2.size(); j++)
			{
				if(j == 0) table.set(j, concat(table.get(j),table2.get(j)));
				else 
				{
					if(idmap.get(table2.get(j)[0]) != null)  //id exists in table (across all tables thus far)
						table.set(j, concat(table.get(idmap.get(table2.get(j)[0])),table2.get(j)));
				}
			}
		}
		StringBuilder sb = new StringBuilder(1000000);
		int good = 0;
		for(String[] str: table)
		{
			
			if(str.length >= cols){
				for(int i = 0; i < str.length-1; i++)
					sb.append(str[i]+"\t");
				sb.append(str[str.length-1]+System.lineSeparator());
				good++;
			}
				
		}
		System.err.println("Final table size: "+good);
		System.out.print(sb);
	}
	
	private String[] concat(String[] s1, String[] s2) {
		String[] comb = new String[s1.length+s2.length-1];
		for(int i = 0; i < s1.length; i++)
			comb[i]=s1[i];
		for(int i = s1.length; i < s1.length+s2.length-1; i++)
			comb[i]=s2[i-s1.length+1];
		return comb;
	}

	private LinkedList<String[]> loadTable(String string) {
		LinkedList<String[]> table = new LinkedList<String[]>();
		File f = new File(string);
		try {
			Scanner s = new Scanner(f);
			while(s.hasNextLine())
			{
				table.add(s.nextLine().split("\t"));
				
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return table;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new tablejoin(args);

	}


}