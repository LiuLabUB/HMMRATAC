package GenomeFileReaders;

/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

import java.awt.Point;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

import Node.TagNode;

public class bedFileReader {
	
	private File file;
	private ArrayList<TagNode> data;
	private int minQ=0;
	private boolean rmDup = false;
	
	
	/**
	 * Constructor which reads data
	 * @param f a BED File containing the data
	 */
	public bedFileReader(File f){
		file = f;
		setData();
	}
	/**
	 * Constructor that reads data. Uses name of file
	 * @param f String that represents the name of the BED File
	 */
	public bedFileReader(String f){
		file = new File(f);
		setData();
	}
	/**
	 * Method for determining which regions of the genome are mappable 
	 * @return an ArrayList of TagNode's representing the mappable regions
	 */
	public ArrayList<TagNode> getMappable(){
		ArrayList<TagNode> mappable = new ArrayList<TagNode>();
		
		Node.TagNode temp = null;
		
		HashMap<String,Point> chrMap = new HashMap<String,Point>();
		for(int i = 0;i < data.size();i++){
			
					
				String chr = data.get(i).getChrom();
				int start = data.get(i).getStart();
				int stop = data.get(i).getStop();
				if (chrMap.containsKey(chr)){
					Point coords = chrMap.get(chr);
					int a = coords.x;
					int b = coords.y;
					if (start < a){
						coords.x = start;
					}
					if (stop > b){
						coords.y = stop;
					}
					chrMap.put(chr, coords);
				}
				else if(!chrMap.containsKey(chr)){
					Point coords = new Point(start,stop);
					chrMap.put(chr, coords);
				}
				
		}
		
		for (String chr : chrMap.keySet()){
			Point p = chrMap.get(chr);
			temp = new TagNode(chr,p.x,p.y);
			mappable.add(temp);
		}	
		return mappable;
	}
	/**
	 * Access the data
	 * @return an ArrayList of TagNode representing the data 
	 */
	public ArrayList<TagNode> getData(){
		return data;
	}
	/**
	 * Access the accepted BED File formats
	 * @return a String that represents the BED File formats that can be parsed
	 */
	public String bedFormats(){
		return "BED3: <chr> <start> <stop>"+"\n"+"BED6: <chr> <start> <stop> <name> <mapQualScore> <strand>"+"\n"+
				"BEDPE: <chr1> <start1> <stop1> <chr2> <start2> <stop2> <name> <mapQ> <strand1> <strand2>";
	}
	/**
	 * Read the input file and set the data
	 */
	private void setData(){
		data = new ArrayList<TagNode>();
		TagNode temp = null;
		Scanner inFile = null;
		try{
			inFile = new Scanner((Readable) new FileReader(file));
		}
		catch(FileNotFoundException e){
			e.printStackTrace();
		}
		int counter = 0;
		while (inFile.hasNextLine()){
			String line = inFile.nextLine();
			String[] fields = line.split("\\s+");
			//System.out.println("length of fields"+"\t"+fields.length);
			if (fields.length == 3 && isBED(fields)){
				temp = readMinBED(fields);
				data.add(temp);
			}
			else if (fields.length > 6 && bedpe(fields) ){
				temp = readBEDPE(fields);
				data.add(temp);
			}
			else if (fields.length == 6 && bedSix(fields)){
				temp = readBEDSix(fields);
				if (counter % 1000000 == 0){System.out.println("reading bedsix");}
				data.add(temp);
			}
			else if (isBED(fields)){
				temp = readMinBED(fields);
				data.add(temp);
			}
			//if (counter % 1000000 == 0){System.out.println(counter);}
			counter+=1;
		}
		if (rmDup == true){
			data = removeDup();
		}
	}
	/**
	 * Method for reading BED-six formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readBEDSix(String[] line){
		TagNode temp = null;
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		String name = line[3];
		int mapQ = Integer.parseInt(line[4]);
		char s1 = line[5].charAt(0);
		if (mapQ >= minQ){
			temp = new TagNode(chr,start,stop,name,mapQ,s1);
		}
		return temp;
	}
	/**
	 * Method for reading BEDPE formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readBEDPE(String[] line){
		TagNode temp = null;
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		int start2 = Integer.parseInt(line[4]);
		int stop2 = Integer.parseInt(line[5]);
		int mapQ = Integer.parseInt(line[7]);
		
		if (mapQ >= minQ){
			temp = new TagNode(chr,Math.min(start, start2),Math.max(stop, stop2));
		}
		return temp;
		
	}
	/**
	 * Method for reading BED-three formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readMinBED(String[] line){
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		TagNode temp = new TagNode(chr,start,stop);
		return temp;
	}
	
	/**
	 * Method to remove duplicate entries in BED file
	 * @return an ArrayList of TagNode representing the data with duplicates removed
	 */
	private ArrayList<TagNode> removeDup(){
		if (rmDup == true){
			HashSet<String> map = new HashSet<String>();
			for (int i = 0; i < data.size();i++){
				String chr = data.get(i).getChrom();
				int start = data.get(i).getStart();
				int stop = data.get(i).getStop();
				String key = chr+"_"+start+"_"+stop;
				map.add(key);
			}
			ArrayList<TagNode> temp = new ArrayList<TagNode>();
			for (String key: map){
				String[] line = key.split("_");
				String chr = line[0];
				int start = Integer.parseInt(line[1]);
				int stop = Integer.parseInt(line[2]);
				TagNode tempNode = new TagNode(chr,start,stop);
				temp.add(tempNode);
			}
			return temp;
		}
		else{return null;}
	}
	
	/**
	 * Determine if the file is in BEDPE format
	 * @param line an Array of Strings representing a line in the BED File
	 * @return a boolean representing whether the data is in BEDPE format
	 */
	private boolean bedpe(String[] line){
		try {
			if (isBED(line)){
				String[] temp = new String[3];
				temp[0] = line[3];temp[1]=line[4];temp[2]=line[5];
				if (isBED(temp) && Math.round(Integer.parseInt(line[7])) == Integer.parseInt(line[7]) 
						&& line[8].equals("+") || line[8].equals("-") 
						&& line[9].equals("+") || line[9].equals("-")){
						return true;
				}
				else{return false;}
			}
			else{return false;}
		}
		catch(NumberFormatException e){
			return false;
		}
	}
	/**
	 * Determine if the file is in BED-six format
	 * @param line an Array of Strings representing a line in the BED File
	 * @return a boolean representing whether the data is in BED-six format
	 */
	private boolean bedSix(String[] line){
		try {
			if (isBED(line) && Math.round(Integer.parseInt(line[4])) == Integer.parseInt(line[4]) && (line[5].equals("+"))
					|| line[5].equals("-")){
				return true;
			}
			else{return false;}
		}
		catch(NumberFormatException e){
			return false;
		}
	}
	
	/**
	 * Determine if the file is in BED format
	 * @param line an Array of Strings representing a line in the BED File
	 * @return a boolean representing whether the data is in BED format
	 */
	private boolean isBED(String[] line){
		try{
			if (line[0].startsWith("chr") || line[0].contains("X") || line[0].contains("Y") || line[0].contains("M")
					|| Math.round(Integer.parseInt(line[0])) == Integer.parseInt(line[0])
					&& Math.round(Integer.parseInt(line[1]))==Integer.parseInt(line[1])
					&& Math.round(Integer.parseInt(line[2]))==Integer.parseInt(line[2])){
				return true;
			
			}
			else{return false;}
		}
		catch (NumberFormatException e){
			return false;
		}
	}
	
	

	
}
