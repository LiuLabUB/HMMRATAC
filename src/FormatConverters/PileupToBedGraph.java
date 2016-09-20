package FormatConverters;

import java.util.ArrayList;
import java.util.HashMap;

import Node.PileupNode2;
import Node.StateNode;
import Node.StateNode2;
import Node.TagNode;

public class PileupToBedGraph {

	
	private ArrayList<TagNode> _bedGraph;
	private ArrayList<PileupNode2> _pileup;
	
	public PileupToBedGraph(ArrayList<PileupNode2> pile){
		_pileup = pile;
		run();
		
	}
	public PileupToBedGraph(ArrayList<PileupNode2> pile,int x){
		_pileup = pile;
		run2();
	}
	
	private void run(){
		HashMap<String,ArrayList<PileupNode2>> map = toMap();
		_bedGraph = new ArrayList<TagNode>();
		for (String chr : map.keySet()){
			_bedGraph.addAll(generateGraph(map.get(chr)));
		}
	}
	
	private ArrayList<TagNode> generateGraph(ArrayList<PileupNode2> states){
		ArrayList<TagNode> bedGraph = new ArrayList<TagNode>();
		TagNode temp = null;
		
		String chr = states.get(0).getChrom();
		int value = states.get(0).getScore2();//needs to change if dealing with a floating score
		int start = states.get(0).getBase();
		int stop;
		double prob = states.get(0).getProb();
		for (int i = 0;i < states.size()-1;i++){
			int newValue = states.get(i).getScore2();//needs to change if dealing with a floating score
			int currentBase = states.get(i).getBase();
			int nextBase = states.get(i+1).getBase();
			if (newValue == value && currentBase==(nextBase-1)){
				if (i > 0){prob += states.get(i).getProb();}
				if (i == states.size()-2){
					stop = states.get(i).getBase();
					temp = new TagNode(chr,start,stop,value);
					temp.setScore3(prob);
					bedGraph.add(temp);
				}
				
			}
			else if (newValue != value || currentBase!=(nextBase-1)){
				stop = states.get(i).getBase();
				temp = new TagNode(chr,start,stop,value);
				temp.setScore3(prob);
				bedGraph.add(temp);
				start = states.get(i).getBase();
				value = newValue;
				prob = states.get(i).getProb();
			}
		}
		return bedGraph;
	}
	
	/*
	 * Generic method for converting pileup to bedgraph
	 */
	private void run2(){
		HashMap<String,ArrayList<PileupNode2>> map = toMap();
		_bedGraph = new ArrayList<TagNode>();
		for (String chr : map.keySet()){
			_bedGraph.addAll(generateGraph2(map.get(chr)));
		}
	}
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
			if (newValue == value && currentBase==(nextBase-1)){
				
				if (i == states.size()-1){
					stop = states.get(i).getBase();
					temp = new TagNode(chr,start,stop,value);
					
					bedGraph.add(temp);
				}
				
			}
			else if (newValue != value || currentBase!=(nextBase-1)){
				stop = states.get(i).getBase();
				temp = new TagNode(chr,start,stop,value);
				
				bedGraph.add(temp);
				start = states.get(i).getBase();
				value = newValue;
				
			}
		}
		return bedGraph;
	}
	
	public ArrayList<TagNode> getBedGraph(){
		return _bedGraph;
	}
	
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
