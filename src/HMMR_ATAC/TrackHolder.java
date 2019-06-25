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
import java.util.List;

import Node.TagNode;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;

public class TrackHolder {
	
	private ArrayList<double[]> tracks;
	private ArrayList<TagNode> positions;
	
	/**
	 * Constructor for creating a TrackHolder object
	 * @param t an ArrayList of doubles representing the data to be held
	 * @param trim an integer representing the number of data points to trim from the left side of the matrix
	 */
	public TrackHolder(ArrayList<double[]> t,int trim){
		tracks = trim(t,trim);
		//positions = pos;
	}
	/**
	 * Access the data as an ArrayList of double arrays
	 * @return An ArrayList of double arrays representing the data
	 */
	public ArrayList<double[]> getRawData(){return tracks;}
	/**
	 * Access the positions
	 */
	public ArrayList<TagNode> getPositions() {return positions;}
	/**
	 * Trim the data
	 * @param t an ArrayList of doubles representing the original data
	 * @param trim an integer representing how many data points to trim
	 * @return an ArrayList of doubles representing the trimmed data
	 */
	public ArrayList<double[]> trim(ArrayList<double[]> t,int trim){
		ArrayList<double[]> updated = new ArrayList<double[]>();
		for (int i = 0;i < t.size();i++){
			double [] temp = t.get(i);
			double [] temp2 = new double[temp.length-trim];
			for (int a = 0; a < temp.length-trim;a++){
				temp2[a] = temp[a];
			}
			updated.add(temp2);
		}
		return updated;
	}
	/**
	 * Access the data as a Dataset, for kmeans
	 * @return a Dataset representing the data for kmeans and javaml applications
	 */
	public Dataset getDataSet(){
		Dataset data = new DefaultDataset();
		for (int i = 0;i < tracks.size();i++){
			DenseInstance ins = new DenseInstance(tracks.get(i));
			//for (int a = 0;a < tracks.get(i).length;a++){
				//System.out.println(tracks.get(i)[a]);
			//}
			data.add(ins);
		}
		
		return data;
	}
	/**
	 * Access the data as a List of List of ObservationVector for baum welch applications
	 * @return a List of List of ObservationVector for baum welch applications 
	 */
	public List<List<ObservationVector>> getBWObs(){
		List<List<ObservationVector>> newList = new ArrayList<List<ObservationVector>>();
		List<ObservationVector> obsList = getObs();
		int a; int i;
		int halfSize = obsList.size()/2;
		
		for (i = 0;i < obsList.size()-1;i+=halfSize){
			
			List<ObservationVector> temp = new ArrayList<ObservationVector>();
			for (a = i;a < i+halfSize-1;a++){
				
				ObservationVector o = (ObservationVector) obsList.get(a);
				temp.add(o);
			}
			newList.add(temp);
		}
		return newList;
	}
	/**
	 * Access the data as a List of ObservationVector for viterbi applications
	 * @return a List of ObservationVector for viterbi applications
	 */
	public List<ObservationVector> getObs(){
		List<ObservationVector> obs = new ArrayList<ObservationVector>();
		for (int i = 0;i < tracks.size();i++){
			ObservationVector vec = new ObservationVector(tracks.get(i));
			obs.add(vec);
		}
		
		return obs;
	}
}
