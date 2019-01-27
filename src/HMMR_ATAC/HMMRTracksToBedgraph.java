package HMMR_ATAC;

import java.util.ArrayList;

import FormatConverters.PileupToBedGraph;
import Node.PileupNode2;
import Node.TagNode;

public class HMMRTracksToBedgraph {
	
	private ArrayList<double[]> tracks;
	private TagNode interval;
	private int step;
	
	private ArrayList<TagNode> nfr;
	private ArrayList<TagNode> mono;
	private ArrayList<TagNode> di;
	private ArrayList<TagNode> tri;
	
	public HMMRTracksToBedgraph(ArrayList<double[]> t, TagNode i,int s){
		tracks = t;
		interval = i;
		step = s;
		run();
	}
	
	public ArrayList<TagNode> getShort(){return nfr;}
	public ArrayList<TagNode> getMono(){return mono;}
	public ArrayList<TagNode> getDi(){return di;}
	public ArrayList<TagNode> getTri(){return tri;}
	
	private void run(){
		ArrayList<TagNode> nfr = new ArrayList<TagNode>();
		ArrayList<TagNode> mono = new ArrayList<TagNode>();
		ArrayList<TagNode> di = new ArrayList<TagNode>();
		ArrayList<TagNode> tri = new ArrayList<TagNode>();
		
		nfr.addAll(runOneCol(0));
		mono.addAll(runOneCol(1));
		di.addAll(runOneCol(2));
		tri.addAll(runOneCol(3));
		
	}

	private ArrayList<TagNode> runOneCol(int c){
		ArrayList<TagNode> temp = new ArrayList<TagNode>();
		int start = interval.getStart();
		ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
		int remainder = interval.getLength() % step;
		int i;
		for ( i = 0;i < tracks.size()-1;i++){
			PileupNode2 pNode = new PileupNode2(start+(i*step),tracks.get(i)[c],interval.getChrom());
			pile.add(pNode);
		}
		PileupNode2 pNode = new PileupNode2(start+(((i)*step)-remainder),tracks.get(i)[c],interval.getChrom());
		pile.add(pNode);
		temp.addAll(new PileupToBedGraph(pile,step).getBedGraph());
		return temp;
	}
}
