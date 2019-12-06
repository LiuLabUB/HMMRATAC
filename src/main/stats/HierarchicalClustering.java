package stats;
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

import org.apache.commons.math3.ml.distance.EuclideanDistance;

public class HierarchicalClustering {
	
	private ArrayList<double[]> input;
	private ArrayList<ClusterNode> clusters;
	private int numClusters = 4;
	
	public HierarchicalClustering(ArrayList<double[]> in){
		input = in;
		createInitial();
		input = null;
	}
	public HierarchicalClustering(ArrayList<double[]> in,int n){
		input = in;
		numClusters = n;
		createInitial();
		input = null;
	}
	
	private void iterate(){
		ArrayList<ClusterNode> temp = new ArrayList<ClusterNode>();
		EuclideanDistance ed = new EuclideanDistance();
		for (int i = 0; i < clusters.size();i++){
			double min = Double.POSITIVE_INFINITY;
			int best = -1;
			for (int a = 0; a < clusters.size();a++){
				if (i != a){
					double dis = ed.compute(clusters.get(i).getKey(), clusters.get(a).getKey());
					if (dis < min){
						min = dis;
						best = a;
					}
				}
			}
			
		}
	}
	
	private void createInitial(){
		clusters = new ArrayList<ClusterNode>();
		for (int i = 0;i < input.size();i++){
			double[] key = input.get(i);
			ArrayList<double[]> values = new ArrayList<double[]>();
			values.add(key);
			clusters.add(new ClusterNode(key,values));
		}
		
	}

}


class ClusterNode{
	private double[] key;
	private ArrayList<double[]> values;
	
	public ClusterNode(double[] k, ArrayList<double[]> v){
		key = k;
		values = v;
	}
	public double[] getKey(){return key;}
	public ArrayList<double[]> getValues(){return values;}
	public void setKey(double[] k){key = k;}
	public void setValues(ArrayList<double[]> v){values = v;}
	public void incrementValues(ArrayList<double[]> v){ values.addAll(v);}
	public void incrementKet(double[] k){
		for (int i = 0;i < key.length;i++){
			key[i] = (key[i] + k[i]) / 2;
		}
	}
}
