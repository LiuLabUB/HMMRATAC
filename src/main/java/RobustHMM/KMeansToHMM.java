package RobustHMM;
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
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.stat.correlation.StorelessCovariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.Instance;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;

public class KMeansToHMM {

	private  Hmm<?> hmm = null;
	/**
	 * Constructor for creating new KMeansToHMM object
	 * @param d a Dataset containing the data
	 * @param K an integer representing the number of states to cluster
	 * @param numIter an integer representing the number of iterations for kmeans clustering
	 * @param diag a boolean to determine whether the resulting covariance matrix should be diagonal
	 * @param equal a boolean to determine whether the initial probability vector should be equal
	 * @param equal2 a boolean to determine whether the transition probability matrix should be equal
	 */
	@SuppressWarnings("unchecked")
	public KMeansToHMM(Dataset d,int K,int numIter,boolean diag,boolean equal,boolean equal2){
		build(d,K,numIter, diag, equal,equal2);
		sort((Hmm<ObservationVector>) hmm);
		
	}
	/**
	 * Sort the hmm
	 * @param hmm an hmm to sort
	 */
	public void sort(Hmm<ObservationVector> hmm) {
		for (int i=0; i<hmm.nbStates(); i++) {
			int indexOfSmallest = indexOfSmallest(hmm,i);
			//E tmp = array[i];
			OpdfMultiGaussian tmp = (OpdfMultiGaussian) hmm.getOpdf(i);
			//array[i] = array[indexOfSmallest];
			hmm.setOpdf(i, hmm.getOpdf(indexOfSmallest));
			hmm.setOpdf(indexOfSmallest, tmp);
			//array[indexOfSmallest] = tmp;
		}
	}
	/**
	 * Find the smallest state, ie the state with the smallest emission values
	 * @param hmm a hmm to find the smallest state
	 * @param startingIndex an integer representing the starting value
	 * @return
	 */
	private int indexOfSmallest(Hmm<ObservationVector> hmm, int startingIndex) {
		int indexOfSmallestSoFar = startingIndex;
		for (int i=startingIndex+1; i<hmm.nbStates(); i++) {
			OpdfMultiGaussian one = (OpdfMultiGaussian) hmm.getOpdf(i);
			OpdfMultiGaussian two = (OpdfMultiGaussian) hmm.getOpdf(indexOfSmallestSoFar);
			if (one.mean()[0] < two.mean()[0]) {
				indexOfSmallestSoFar = i;
			}
		}
		return indexOfSmallestSoFar;
	}
	/**
	 * Access the completed hmm
	 * @return the completed hmm
	 */
	public Hmm<?> getHMM(){return hmm;}
	/**
	 * Build the model with kmeans clustering
	 * @param data a Dataset containing the data
	 * @param K an integer representing the number of states to cluster
	 * @param numIter an integer representing the number of iterations for kmeans clustering
	 * @param diag a boolean to determine whether the resulting covariance matrix should be diagonal
	 * @param equal a boolean to determine whether the initial probability vector should be equal
	 * @param equal2 a boolean to determine whether the transition probability matrix should be equal
	 */
	private void build(Dataset data,int K,int numIter,boolean diag,boolean equal,boolean equal2){
		int numFeatures = data.noAttributes();
		
		KMeans k = new KMeans(K,numIter);//changed to Kmeans (modified version with uniform centroids)
		//k.setUniformInitialCentroids();//added this line when changed to modified version
		Dataset[] clustered = k.cluster(data);
		int[] assignments = new int[data.size()];
		List<OpdfMultiGaussian> opdf = new ArrayList<OpdfMultiGaussian>();
		
		for (int a = 0; a < clustered.length;a++){
			Dataset cluster = clustered[a];
			double[] means = new double[numFeatures];
			StorelessCovariance cov = new StorelessCovariance(numFeatures);
			double[][] var = new double[numFeatures][cluster.size()];
			for (int x = 0;x < cluster.size();x++){
				Instance ins = cluster.get(x);
				assignments[ins.getID()] = a;
				Iterator<Double> iter = ins.values().iterator();
				double[] values = new double[numFeatures];
				int counter = 0;
				while(iter.hasNext()){
					double value = iter.next();
					values[counter] = value;
					means[counter] += value;
					var[counter][x] = value;
					counter++;
				}
				cov.increment(values);
				
			}
			for (int y = 0;y < means.length;y++){
				means[y] /= (double) cluster.size();
			}
			double[][] covMat = null;
			if (diag == false){
				covMat = cov.getData();
			}
			else{
				covMat = new double[numFeatures][numFeatures];
				Variance variance = new Variance();
				for (int z = 0;z<numFeatures;z++){
					covMat[z][z] = variance.evaluate(var[z]);
				}
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(means, covMat);
			opdf.add(pdf);
		}
		
		double[][] trans = new double[K][K];
		if (!equal2){
			for (int a = 0;a < assignments.length-1;a++){
				trans[assignments[a]][assignments[a+1]]++;
			}
			for (int i = 0; i < trans.length;i++){
			
				int elementCounter=0;
				for (int x = 0;x < trans[i].length;x++){
					elementCounter += trans[i][x];
				}
				for (int y = 0; y < trans[i].length;y++){
					trans[i][y] /= elementCounter;
				}
			}
		}
		else{
			for (int i = 0;i < trans.length;i++){
				for (int a = 0;a < trans[i].length;a++){
					trans[i][a] = 1.0/(double)K;
				}
			}
		}
		
		double[] initial = new double[K];
		if (equal == true){
			for (int i = 0; i < K;i++){
				initial[i] = 1/(double)K;
			}
		}
		hmm = new Hmm<ObservationVector>(initial,trans,opdf);
	}
	
	
	
	
	
}
