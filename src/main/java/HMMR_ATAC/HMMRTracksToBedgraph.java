package HMMR_ATAC;

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
