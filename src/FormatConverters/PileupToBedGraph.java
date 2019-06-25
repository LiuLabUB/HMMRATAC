package FormatConverters;

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

import java.util.ArrayList;
import java.util.HashMap;

import Node.PileupNode2;
import Node.TagNode;

public class PileupToBedGraph {

	
	private ArrayList<TagNode> _bedGraph;
	private ArrayList<PileupNode2> _pileup;
	private int step = 1;
	
	/**
	 * Constructor for creating new PileupToBedGraph object and generating the data
	 * @param pile an ArrayList of PileupNode2's to convert into a bedgraph
	 * @param x an integer representing the step in the pileup
	 */
	public PileupToBedGraph(ArrayList<PileupNode2> pile,int x){
		_pileup = pile;
		step = x;
		run();
	}
	
	
	
	/**
	 * Convert the entire pileup into bedgraph
	 */
	private void run(){
		HashMap<String,ArrayList<PileupNode2>> map = toMap();
		_bedGraph = new ArrayList<TagNode>();
		for (String chr : map.keySet()){
			_bedGraph.addAll(generateGraph2(map.get(chr)));
		}
	}
	/**
	 * Convert one portion of the pileup into bedgraph
	 * @param states an ArrayList of PileupNode2's to convert
	 * @return an ArrayList of TagNode's representing the bedgraph data
	 */
	private ArrayList<TagNode> generateGraph2(ArrayList<PileupNode2> states){
		ArrayList<TagNode> bedGraph = new ArrayList<TagNode>();
		TagNode temp = null;
		
		String chr = states.get(0).getChrom();
		double value = states.get(0).getScore();
		int start = states.get(0).getBase();
		int stop;
		
		for (int i = 0;i < states.size()-1;i++){
			double newValue = states.get(i).getScore();
			int currentBase = states.get(i).getBase();
			int nextBase = states.get(i+1).getBase();
			if (newValue == value && currentBase==(nextBase-step)){
				
				if (i == states.size()-1){
					stop = states.get(i).getBase();
					temp = new TagNode(chr,start,stop,value);
					
					bedGraph.add(temp);
				}
				
			}
			else if (newValue != value || currentBase!=(nextBase-step)){
				stop = states.get(i).getBase();
				temp = new TagNode(chr,start,stop,value);
				
				bedGraph.add(temp);
				start = states.get(i).getBase();
				value = newValue;
				
			}
		}
		return bedGraph;
	}
	/**
	 * Access the finished bedgraph data
	 * @return an ArrayList of TagNode's representing the completed bedgraph
	 */
	public ArrayList<TagNode> getBedGraph(){
		return _bedGraph;
	}
	/**
	 * Split the pileup into individual chromosomes for ease of computation
	 * @return a HashMap of ArrayList's of PileupNode2's. Each entry in the HasHMap represents one chromosome
	 */
	private HashMap<String,ArrayList<PileupNode2>> toMap(){
		HashMap<String,ArrayList<PileupNode2>> map = new HashMap<String,ArrayList<PileupNode2>>();
		
		for (int i = 0; i < _pileup.size();i++){
			String key = _pileup.get(i).getChrom();
			if (map.containsKey(key)){
				ArrayList<PileupNode2> temp = map.get(key);
				temp.add(_pileup.get(i));
				map.put(key, temp);
			}
			else if (!map.containsKey(key)){
				ArrayList<PileupNode2> temp = new ArrayList<PileupNode2>();
				temp.add(_pileup.get(i));
				map.put(key, temp);
			}
		}
		
		return map;
	}
}
