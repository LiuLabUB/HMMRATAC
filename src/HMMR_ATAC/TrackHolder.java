package HMMR_ATAC;

import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;

public class TrackHolder {
	
	private ArrayList<double[]> tracks;
	
	/**
	 * Constructor for creating a TrackHolder object
	 * @param t an ArrayList of doubles representing the data to be held
	 * @param trim an integer representing the number of data points to trim from the left side of the matrix
	 */
	public TrackHolder(ArrayList<double[]> t,int trim){
		tracks = trim(t,trim);
	}
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
