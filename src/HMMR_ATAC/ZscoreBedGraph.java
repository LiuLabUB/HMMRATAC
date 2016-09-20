package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class ZscoreBedGraph {
	
	private ArrayList<TagNode> bedgraph;
	private double mean;
	private double std;
	private double cutoff;
	
	private ArrayList<TagNode> output;
	
	public ZscoreBedGraph(ArrayList<TagNode> b,double m,double s,double c){
		bedgraph = b;
		mean = m;
		std = s;
		cutoff = c;
		locate();
	}
	public ArrayList<TagNode> getResults(){return output;}
	private void locate(){
		output = new ArrayList<TagNode>();
		for (int i = 0; i < bedgraph.size();i++){
			double value = bedgraph.get(i).getScore2();
			double z = (value - mean)/std;
			if (z >= cutoff){
				output.add(bedgraph.get(i));
			}
		}
	}

}
