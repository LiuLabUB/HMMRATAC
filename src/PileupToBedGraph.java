package FormatConverters;

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
