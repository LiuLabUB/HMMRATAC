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
