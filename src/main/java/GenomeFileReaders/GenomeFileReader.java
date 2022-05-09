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

	private final ArrayList<TagNode> _genMap = new ArrayList<>();
	/**
	 * Constructor for creating Object and reading the data
	 * @param input a File representing the genome data
	 */
	public GenomeFileReader(File input) throws FileNotFoundException{
		Scanner inFile =new Scanner(new FileReader(input));
		
		while (inFile.hasNext()) {
			String chr = inFile.next();
			int size = inFile.nextInt();
			_genMap.add(new TagNode(chr,0,size));
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
