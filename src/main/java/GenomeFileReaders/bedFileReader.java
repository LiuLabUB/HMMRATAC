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
	
	private final File file;
	private ArrayList<TagNode> data;
	private final int minQ = 0;
	private boolean rmDup = false;

	/**
	 * Constructor which reads data
	 * @param file a BED File containing the data
	 * @param rmDup boolean flag whether to remove duplicate entries in the BED file
	 */
	public bedFileReader(File file, boolean rmDup) {
		this.rmDup = rmDup;
		this.file = file;
		setData();
	}

	/**
	 * Constructor which reads data
	 * @param file a BED File containing the data
	 */
	public bedFileReader(File file) {
		this.file = file;
		setData();
	}
	/**
	 * Constructor that reads data. Uses name of file
	 * @param f String that represents the name of the BED File
	 */
	public bedFileReader(String f) {
		this(new File(f));
	}
	/**
	 * Method for determining which regions of the genome are mappable 
	 * @return an ArrayList of TagNode's representing the mappable regions
	 */
	public ArrayList<TagNode> getMappable(){
		ArrayList<TagNode> mappable = new ArrayList<>();
		HashMap<String,Point> chrMap = new HashMap<>();

		for (TagNode datum : data) {

			String chr = datum.getChrom();
			int start = datum.getStart();
			int stop = datum.getStop();
			Point coords;
			if (chrMap.containsKey(chr)) {
				coords = chrMap.get(chr);
				coords.x = Math.min(start, coords.x);
				coords.y = Math.max(stop, coords.y);
			} else {
				coords = new Point(start, stop);
			}
			chrMap.put(chr, coords);
		}
		
		for (String chr : chrMap.keySet()){
			Point p = chrMap.get(chr);
			mappable.add(new TagNode(chr,p.x,p.y));
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
	public String bedFormats() {
		return "BED3: <chr> <start> <stop>\nBED6: <chr> <start> <stop> <name> <mapQualScore> <strand>\n" +
				"BEDPE: <chr1> <start1> <stop1> <chr2> <start2> <stop2> <name> <mapQ> <strand1> <strand2>";
	}

	/**
	 * Read the input file and set the data
	 */
	private void setData() {
		data = new ArrayList<>();
		try {
			Scanner inFile = new Scanner(new FileReader(file));
			while (inFile.hasNextLine()) {
				String line = inFile.nextLine();
				String[] fields = line.split("\\s+");

				if (fields.length == 3 && isBED(fields)) {
					data.add(readMinBED(fields));
				} else if (fields.length > 6 && bedpe(fields)) {
					data.add(readBEDPE(fields));
				} else if (fields.length == 6 && bedSix(fields)) {
					data.add(readBEDSix(fields));
				} else if (isBED(fields)) {
					data.add(readMinBED(fields));
				}
			}
			if (rmDup) {
				data = removeDup();
			}
		}
		catch(FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Method for reading BED-six formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readBEDSix(String[] line) {
		TagNode node = null;
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		String name = line[3];
		int mapQ = Integer.parseInt(line[4]);
		char s1 = line[5].charAt(0);
		if (mapQ >= minQ){
			node = new TagNode(chr, start, stop, name, mapQ, s1);
		}
		return node;
	}

	/**
	 * Method for reading BEDPE formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readBEDPE(String[] line) {
		TagNode node = null;
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		int start2 = Integer.parseInt(line[4]);
		int stop2 = Integer.parseInt(line[5]);
		int mapQ = Integer.parseInt(line[7]);
		
		if (mapQ >= minQ) {
			node = new TagNode(chr, Math.min(start, start2), Math.max(stop, stop2));
		}
		return node;
		
	}

	/**
	 * Method for reading BED-three formatted files
	 * @param line an Array of Strings representing a line in the BED File
	 * @return	a TagNode representing the BED data
	 */
	private TagNode readMinBED(String[] line) {
		String chr = line[0];
		int start = Integer.parseInt(line[1]);
		int stop = Integer.parseInt(line[2]);
		return new TagNode(chr,start,stop);
	}
	
	/**
	 * Method to remove duplicate entries in BED file
	 * @return an ArrayList of TagNode representing the data with duplicates removed
	 */
	private ArrayList<TagNode> removeDup(){
		HashSet<String> map = new HashSet<>();
		for (TagNode datum : data) {
			String chr = datum.getChrom();
			int start = datum.getStart();
			int stop = datum.getStop();
			String key = chr + "_" + start + "_" + stop;
			map.add(key);
		}
		ArrayList<TagNode> list = new ArrayList<>();
		for (String key: map){
			String[] line = key.split("_");
			String chr = line[0];
			int start = Integer.parseInt(line[1]);
			int stop = Integer.parseInt(line[2]);
			TagNode tempNode = new TagNode(chr, start, stop);
			list.add(tempNode);
		}
		return list;
	}
	
	/**
	 * Determine if the file is in BEDPE format
	 * @param line an Array of Strings representing a line in the BED File
	 * @return a boolean representing whether the data is in BEDPE format
	 */
	private boolean bedpe(String[] line) {
		try {
			if (isBED(line)) {
				String[] arr = {line[3], line[4], line[5]};

				return isBED(arr) && Math.round(Integer.parseInt(line[7])) == Integer.parseInt(line[7])
						&& line[8].equals("+") || line[8].equals("-")
						&& line[9].equals("+") || line[9].equals("-");
			}
			else {return false;}
		}
		catch (NumberFormatException e) {
			return false;
		}
	}

	/**
	 * Determine if the file is in BED-six format
	 * @param line an Array of Strings representing a line in the BED File
	 * @return a boolean representing whether the data is in BED-six format
	 */
	private boolean bedSix(String[] line) {
		try {
			return isBED(line) && Math.round(Integer.parseInt(line[4])) == Integer.parseInt(line[4])
					&& (line[5].equals("+")) || line[5].equals("-");
		}
		catch (NumberFormatException e) {
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
			return line[0].startsWith("chr") || line[0].contains("X") || line[0].contains("Y") || line[0].contains("M")
					|| Math.round(Integer.parseInt(line[0])) == Integer.parseInt(line[0])
					&& Math.round(Integer.parseInt(line[1])) == Integer.parseInt(line[1])
					&& Math.round(Integer.parseInt(line[2])) == Integer.parseInt(line[2]);
		}
		catch (NumberFormatException e) {
			return false;
		}
	}
}
