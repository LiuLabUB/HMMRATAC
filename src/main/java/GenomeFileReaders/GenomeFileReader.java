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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import Node.TagNode;

public class GenomeFileReader {

	private ArrayList<TagNode> _genMap;
	/**
	 * Constructor for creating Object and reading the data
	 * @param input a File representing the genome data
	 * @throws FileNotFoundException
	 */
	public GenomeFileReader(File input) throws FileNotFoundException{
		_genMap = new ArrayList<TagNode>();
		TagNode temp = null;
		@SuppressWarnings("resource")
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		
		while(inFile.hasNext()){
			String chr = inFile.next();
			int size = inFile.nextInt();
			temp = new TagNode(chr,0,size);//Note the 20000 value is so it can be used in Bin.java
			_genMap.add(temp);
			
			
		}
		
	}
	/**
	 * Constructor for creating Object and reading the data
	 * @param input a String representing the name of the file containing the genome data
	 * @throws FileNotFoundException
	 */
	public GenomeFileReader(String input) throws FileNotFoundException{
		_genMap = new ArrayList<TagNode>();
		TagNode temp = null;
		@SuppressWarnings("resource")
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		
		while(inFile.hasNext()){
			String chr = inFile.next();
			int size = inFile.nextInt();
			temp = new TagNode(chr,0,size);//Note the 20000 value is so it can be used in Bin.java
			_genMap.add(temp);
			
			
		}
		
	}
	/**
	 * Access the data
	 * @return an ArrayList of TagNode representing the data
	 */
	public ArrayList<TagNode> getMap(){
		return _genMap;
	}
	
}
