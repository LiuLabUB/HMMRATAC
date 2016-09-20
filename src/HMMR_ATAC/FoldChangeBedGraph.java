package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class FoldChangeBedGraph {
	
	private ArrayList<TagNode> bedgraph;
	private double mean;
	private double upper;
	private double lower;
	
	private ArrayList<TagNode> output;
	
	public FoldChangeBedGraph(ArrayList<TagNode> b, double m, double u,double l){
		bedgraph = b;
		mean = m;
		upper = u;
		lower = l;
		locate();
	}
	public ArrayList<TagNode> getResults(){return output;}
	private void locate(){
		output = new ArrayList<TagNode>();
		for (int i = 0; i < bedgraph.size();i++){
			double fc = bedgraph.get(i).getScore2() / mean;
			if (fc >= lower && fc <= upper){
				output.add(bedgraph.get(i));
			}
		}
	}

}
