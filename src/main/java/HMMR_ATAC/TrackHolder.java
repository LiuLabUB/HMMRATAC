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

import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;

import java.util.ArrayList;
import java.util.List;

public class TrackHolder {
	
	private final ArrayList<double[]> tracks;
	
	/**
	 * Constructor for creating a TrackHolder object
	 *
	 * @param data    an ArrayList of doubles representing the data to be held
	 * @param trim an integer representing the number of data points to trim from the left side of the matrix
	 */
	public TrackHolder(ArrayList<double[]> data, int trim) {
		tracks = trim(data, trim);
	}
	
	/**
	 * Access the data as an ArrayList of double arrays
	 *
	 * @return An ArrayList of double arrays representing the data
	 */
	public ArrayList<double[]> getRawData() {
		return tracks;
	}
	
	/**
	 * Trim the data
	 *
	 * @param data    an ArrayList of doubles representing the original data
	 * @param trim an integer representing how many data points to trim
	 * @return an ArrayList of doubles representing the trimmed data
	 */
	public ArrayList<double[]> trim(ArrayList<double[]> data, int trim) {
		ArrayList<double[]> updated = new ArrayList<>();
		for (double[] temp : data) {
			double[] temp2 = new double[temp.length - trim];
			if (temp.length - trim >= 0) System.arraycopy(temp, 0, temp2, 0, temp.length - trim);
			updated.add(temp2);
		}
		return updated;
	}
	
	/**
	 * Access the data as a Dataset, for kmeans
	 *
	 * @return a Dataset representing the data for kmeans and javaml applications
	 */
	public Dataset getDataSet() {
		Dataset data = new DefaultDataset();
		for (double[] track : tracks) {
			data.add(new DenseInstance(track));
		}
		return data;
	}
	
	/**
	 * Access the data as a List of Lists of ObservationVector for baum welch applications
	 *
	 * @return a List of Lists of ObservationVector for baum welch applications
	 */
	public List<List<ObservationVector>> getBWObs() {
		List<List<ObservationVector>> newList = new ArrayList<>();
		List<ObservationVector> obsList = getObs();
		int halfSize = obsList.size() / 2;
		
		for (int i = 0; i < obsList.size() - 1; i += halfSize) {
			
			List<ObservationVector> temp = new ArrayList<>();
			for (int a = i; a < i + halfSize - 1; a++) {
				
				ObservationVector o = obsList.get(a);
				temp.add(o);
			}
			newList.add(temp);
		}
		return newList;
	}
	
	/**
	 * Access the data as a List of ObservationVector for viterbi applications
	 *
	 * @return a List of ObservationVector for viterbi applications
	 */
	public List<ObservationVector> getObs() {
		List<ObservationVector> obs = new ArrayList<>();
		for (double[] track : tracks) {
			obs.add(new ObservationVector(track));
		}
		return obs;
	}
}
