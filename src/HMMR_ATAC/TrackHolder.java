package HMMR_ATAC;

import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;

public class TrackHolder {
	
	private ArrayList<double[]> tracks;
	
	public TrackHolder(ArrayList<double[]> t){
		tracks = t;
	}

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
	public List<ObservationVector> getObs(){
		List<ObservationVector> obs = new ArrayList<ObservationVector>();
		for (int i = 0;i < tracks.size();i++){
			ObservationVector vec = new ObservationVector(tracks.get(i));
			obs.add(vec);
		}
		
		return obs;
	}
}
