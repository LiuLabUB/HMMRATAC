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
import java.util.List;

import JAHMMTest.BaumWelchScaledLearner;
import JAHMMTest.FitRobust;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;

public class BaumWelch {
	
	private Hmm<?> h;
	private List<List<ObservationVector>> obs;
	private int maxIter;
	private double epsilon;
	/**
	 * Constructor for creating new BaumWelch object
	 * @param H a hidden markov model
	 * @param o a List of List's of ObservationVector's
	 */
	public BaumWelch(Hmm<?> H, List<List<ObservationVector>> o){
		h = H;
		obs = o;
		maxIter=150;
		epsilon = 0.001;
	}
	/**
	 * Constructor for creating new BaumWelch object
	 * @param H a hidden markov model
	 * @param o a List of List of ObservationVector
	 * @param i an integer representing the maximum iterations to perform
	 */
	public BaumWelch(Hmm<?> H, List<List<ObservationVector>> o,int i){
		h=H;
		obs=o;
		maxIter=i;
		epsilon = 0.001;
	}
	/**
	 * Constructor for creating new BaumWelch object
	 * @param H a hidden markov model
	 * @param o a List of List of ObservationVector
	 * @param i an integer representing the maximum iterations to perform
	 * @param e a double representing the epsilon to check for model covergence
	 */
	public BaumWelch(Hmm<?> H, List<List<ObservationVector>> o,int i,double e){
		h=H;
		obs=o;
		maxIter=i;
		epsilon = e;
	}
	
	/**
	 * Build the model
	 * @return a refined model after Baum Welch training
	 */
	@SuppressWarnings("unchecked")
	public Hmm<ObservationVector> build(){
		Hmm<ObservationVector> firstHmm = (Hmm<ObservationVector>) h;
		firstHmm = checkModel(firstHmm);
		BaumWelchScaledLearner sbw = new BaumWelchScaledLearner();
		Hmm<ObservationVector> scaled = null;
		int iter = 0;
		while (iter < maxIter  ){
			scaled = sbw.iterate(firstHmm, obs);
			scaled = checkModel(scaled);
			if (converged(scaled,firstHmm)){
				break;
			}
			iter += 1;
			firstHmm = scaled;
		}
		//Set proportional initial probabilities
		for (int i = 0; i < scaled.nbStates();i++){
			scaled.setPi(i,(double) 1/scaled.nbStates()); //bug. originally hardcoded at 0.25, but should be flexable to other K's
		}
		if (!Double.isNaN(scaled.getAij(0, 0))){
			return scaled;
		}
		else{
			return (Hmm<ObservationVector>) h;
		}
	}
	/**
	 * Check the model 
	 * @param hmm HMM to check
	 * @return modified HMM after robust correction
	 */
	private Hmm<ObservationVector> checkModel(Hmm<ObservationVector> hmm){
		for (int i = 0;i < hmm.nbStates();i++){
			OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
			double[][] cov = pdf.covariance();
			FitRobust fitter = new FitRobust(cov);
			double[][] temp = fitter.getCovariance();
			OpdfMultiGaussian t = new OpdfMultiGaussian(pdf.mean(),temp);
			hmm.setOpdf(i, t);
		}
		return hmm;
	}
	/**
	 * Access whether the model converged
	 * @param h1 HMM before BW iteration
	 * @param h2 HMM after BW iteration
	 * @return	a boolean representing whether the model converged within epsilon value
	 */
	private boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2){
		int counter = 0;
		for (int i = 0;i < h1.nbStates();i++){
			OpdfMultiGaussian pdf1 = (OpdfMultiGaussian) h1.getOpdf(i);
			OpdfMultiGaussian pdf2 = (OpdfMultiGaussian) h2.getOpdf(i);
			double[] value1 = pdf1.mean();
			double[] value2 = pdf2.mean();
			for (int a = 0; a < value1.length;a++){
				if (Math.abs(value1[a] - value2[a]) > epsilon){
					
					counter += 1;
				}
			}
		}
		
		if (counter == 0){
			return true;
		}
		return false;
	}
}
