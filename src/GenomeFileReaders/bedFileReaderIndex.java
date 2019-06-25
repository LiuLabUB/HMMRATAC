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
import java.util.ArrayList;
import java.util.HashMap;

import Node.TagNode;

public class bedFileReaderIndex {
	
	private File input;
	private File genome;
	
	public bedFileReaderIndex(File i, File g){
		input = i;
		genome = g;
	}
	
	private void read() throws FileNotFoundException{
		GenomeFileReader gReader = new GenomeFileReader(genome);
		
		ArrayList<TagNode> gen = gReader.getMap();
		HashMap<String,ArrayList<TagNode>> IndexHash = new HashMap<String,ArrayList<TagNode>>();
		HashMap<Integer,ArrayList<TagNode>> readHash = new HashMap<Integer,ArrayList<TagNode>>();
		
		int counter = 0;
		for (int i = 0;i < gen.size();i++){
			String chr = gen.get(i).getChrom();
			int stop = gen.get(i).getStop();
			ArrayList<TagNode> tempNode = new ArrayList<TagNode>();
			
			for (int a = 0;a < stop;a += 500000){
				TagNode t = new TagNode(chr,a,a+500000,counter);
				tempNode.add(t);
				counter++;
			}
			IndexHash.put(chr, tempNode);
		}
		
		
		
	}

}
