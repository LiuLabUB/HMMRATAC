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

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;

import java.util.ArrayList;
import java.util.List;

public class RobustHMM {
	//inputs
	private final List<List<?>> train;
	private final List<?> obs;
	private final Hmm<?> hmm;
	private final boolean welch;
	private final int k;
	private final String method;
	private final int numDist;
	private int BWIter = 10;
	//specifics
	private KMeansLearner<?> kmeans;
	private final BaumWelchLearner bw;
	private final BaumWelchScaledLearner sbw;
	
	private int[] states;
	private Hmm<?> outHmm;
	private double prob;
	private double prob2;
	
	/**
	 * Constructor
	 *
	 * @param o a List of Observations representing data to be decoded
	 * @param t a List of List of Observation representing training data
	 * @param h an initial hmm
	 * @param w a boolean to determine whether to perform baum welch training
	 * @param K an integer representing the number of states in the model
	 * @param m a String representing the type of data
	 * @param n an integer representing the number of distributions in the Observation data
	 */
	public RobustHMM(List<?> o, List<List<?>> t, Hmm<?> h, boolean w, int K, String m, int n) {
		train = t;
		obs = o;
		hmm = h;
		welch = w;
		k = K;
		method = m;
		numDist = n;
		bw = new BaumWelchLearner();
		sbw = new BaumWelchScaledLearner();
		run();
	}
	
	/**
	 * Set the number of baum welch iterations
	 *
	 * @param n
	 */
	public void setBaumWelchMaxIter(int n) {
		BWIter = n;
	}
	
	/**
	 * Access the state annotations
	 *
	 * @return an Array of integers representing the state assignments after decoding
	 */
	public int[] getStates() {
		return states;
	}
	
	/**
	 * Access state probability
	 *
	 * @return a double representing the state probability
	 */
	public double getProb() {
		return prob;
	}
	
	/**
	 * Access state probability
	 *
	 * @return a double representing the state probability
	 */
	public double getProb2() {
		return prob2;
	}
	
	/**
	 * Access the final hmm
	 *
	 * @return the final hmm
	 */
	public Hmm<?> getHmm() {
		return outHmm;
	}
	
	/**
	 * Run the data
	 */
	private void run() {
		if (isInt()) {
			runInt();
		} else if (isReal()) {
			runReal();
		} else if (isVector()) {
			runVector();
		} else if (isMixture()) {
			runMixture();
		}
	}
	
	/**
	 * Run a mixture style dataset
	 */
	@SuppressWarnings("unchecked")
	private void runMixture() {
		List<List<ObservationReal>> newTrain = formatMixtureTrain();
		List<ObservationReal> newObs = formatRealObs();
		Hmm<ObservationReal> newHmm = null;
		if (k > 0) {
			newHmm = (Hmm<ObservationReal>) kmeans.learn();
		} else {
			newHmm = formatRealHmm();
		}
		if (welch) {
			for (int x = 0; x < BWIter; x++) {
				newHmm = bw.iterate(newHmm, newTrain);
				if (bwConverge(newHmm)) {
					break;
				}
			}
			assert newHmm != null;
			if (!bwConverge(newHmm)) {
				for (int y = 0; y < BWIter; y++) {
					newHmm = sbw.iterate(newHmm, newTrain);
					if (bwConverge(newHmm)) {
						break;
					}
				}
			}
			if (!bwConverge(newHmm)) {
				System.out.println("Baum-Welch Iterations" + "\t" + BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
		}
		states = viterbiReal(newHmm, newObs);
		outHmm = newHmm;
	}
	
	/**
	 * Run a multivariate gaussian style dataset
	 */
	@SuppressWarnings("unchecked")
	private void runVector() {
		List<List<ObservationVector>> newTrain = formatVectorTrain();
		List<ObservationVector> newObs = formatVectorObs();
		Hmm<ObservationVector> newHmm = null;
		if (k > 0) {
			newHmm = (Hmm<ObservationVector>) kmeans.learn();
		} else {
			newHmm = formatVectorHmm();
		}
		if (welch) {
			
			for (int x = 0; x < BWIter; x++) {
				newHmm = sbw.iterate(newHmm, newTrain);
				if (bwConverge(newHmm)) {
					System.out.println("Baum Welch Error!!!");
					System.exit(1);
					break;
				}
			}
		}
		states = viterbiVector(newHmm, newObs);
		outHmm = newHmm;
	}
	
	/**
	 * Run a double style dataset
	 */
	@SuppressWarnings("unchecked")
	private void runReal() {
		List<List<ObservationReal>> newTrain = formatRealTrain();
		List<ObservationReal> newObs = formatRealObs();
		Hmm<ObservationReal> newHmm = null;
		if (k > 0) {
			newHmm = (Hmm<ObservationReal>) kmeans.learn();
		} else {
			newHmm = formatRealHmm();
		}
		if (welch) {
			for (int x = 0; x < BWIter; x++) {
				newHmm = bw.iterate(newHmm, newTrain);
				if (bwConverge(newHmm)) {
					break;
				}
			}
			assert newHmm != null;
			if (!bwConverge(newHmm)) {
				for (int y = 0; y < BWIter; y++) {
					newHmm = sbw.iterate(newHmm, newTrain);
					if (bwConverge(newHmm)) {
						break;
					}
				}
			}
			if (!bwConverge(newHmm)) {
				System.out.println("Baum-Welch Iterations" + "\t" + BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
		}
		states = viterbiReal(newHmm, newObs);
		outHmm = newHmm;
	}
	
	/**
	 * Run a integer style dataset
	 */
	@SuppressWarnings("unchecked")
	private void runInt() {
		List<List<ObservationInteger>> newTrain = formatIntTrain();
		List<ObservationInteger> newObs = formatIntObs();
		Hmm<ObservationInteger> newHmm = null;
		if (k > 0) {
			newHmm = (Hmm<ObservationInteger>) kmeans.learn();
		} else {
			newHmm = formatIntHmm();
		}
		if (welch) {
			for (int x = 0; x < BWIter; x++) {
				newHmm = bw.iterate(newHmm, newTrain);
				if (bwConverge(newHmm)) {
					break;
				}
			}
			assert newHmm != null;
			if (!bwConverge(newHmm)) {
				for (int y = 0; y < BWIter; y++) {
					newHmm = sbw.iterate(newHmm, newTrain);
					if (bwConverge(newHmm)) {
						break;
					}
				}
			}
			if (!bwConverge(newHmm)) {
				System.out.println("Baum-Welch Iterations" + "\t" + BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
		}
		states = viterbiInt(newHmm, newObs);
		outHmm = newHmm;
	}
	
	/**
	 * Run viterbi on integer data
	 *
	 * @param h a hmm
	 * @param l a List of ObservationInteger representing the data
	 * @return an Array of integers representing the state assignments
	 */
	private int[] viterbiInt(Hmm<ObservationInteger> h, List<ObservationInteger> l) {
		ViterbiCalculator vit = new ViterbiCalculator(l, h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	
	/**
	 * Run viterbi on double data
	 *
	 * @param h a hmm
	 * @param l a List of ObservationReal representing the data
	 * @return an Array of integers representing the state assignments
	 */
	private int[] viterbiReal(Hmm<ObservationReal> h, List<ObservationReal> l) {
		ViterbiCalculator vit = new ViterbiCalculator(l, h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	
	/**
	 * Run viterbi on multivariate gaussian data
	 *
	 * @param h a hmm
	 * @param l a List of ObservationVector representing the data
	 * @return an Array of integers representing the state assignments
	 */
	private int[] viterbiVector(Hmm<ObservationVector> h, List<ObservationVector> l) {
		ViterbiCalculator vit = new ViterbiCalculator(l, h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	
	/**
	 * Format the data as multivariate gaussian data
	 *
	 * @return a List of List of ObservationVector representing the data for baum welch training
	 */
	@SuppressWarnings("unchecked")
	private List<List<ObservationVector>> formatVectorTrain() {
		List<List<ObservationVector>> newList = new ArrayList<>();
		if (train != null) {
			for (List<?> objects : train) {
				ArrayList<ObservationVector> tempList = (ArrayList<ObservationVector>) objects;
				newList.add(tempList);
			}
			int dim = newList.get(0).get(0).dimension();
			if (k > 0) {
				kmeans = new KMeansLearner<>(k, new OpdfMultiGaussianFactory(dim), newList);
			}
		}
		/*
		 * This is a test. If no training set is given, duplicate the observation set and use that.
		 * Could also attempt to split up the observation set into smaller, random pieces.
		 * 	Duplicating observation set resulted in out of memory error
		 * 	Splitting the data into random chunks resulted in non convergence of BW
		 * Try splitting the data into 1000bp chunks (for whole data)
		 * Note: either of these approaches may not work.
		 * Update 8/27/15:
		 * 	Splitting into 1000bp chunks sometimes worked. Problem was that although individual tracks were not
		 * 	correlated in entire dataset, they were sometimes within the chunks, leading to problems
		 * 	Try splitting entire dataset in half instead. Less likely for data tracks to be correlated
		 * Original approach was simply to return null if train == null
		 */
		else {
			//split the data into 1000bp chunks
			int a;
			int i;
			int halfSize = obs.size() / 2;
			for (i = 0; i < obs.size() - 1; i += halfSize) {
				List<ObservationVector> temp = new ArrayList<>();
				for (a = i; a < i + halfSize - 1; a++) {
					ObservationVector o = (ObservationVector) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
			
			
		}
		return newList;
	}
	
	/**
	 * Format the data as multivariate gaussian data
	 *
	 * @return a List of ObservationVector representing the data for viterbi decoding
	 */
	@SuppressWarnings("unchecked")
	private List<ObservationVector> formatVectorObs() {
		return (List<ObservationVector>) (ArrayList<ObservationVector>) obs;
	}
	
	/**
	 * Format hmm for multicariate gaussian data
	 *
	 * @return a hmm
	 */
	@SuppressWarnings("unchecked")
	private Hmm<ObservationVector> formatVectorHmm() {
		if (hmm != null) {
			return (Hmm<ObservationVector>) hmm;
		} else {
			return null;
		}
	}
	
	/**
	 * Format the data as mixture data
	 *
	 * @return a List of List of ObservationReal representing the data for baum welch training
	 */
	@SuppressWarnings("unchecked")
	private List<List<ObservationReal>> formatMixtureTrain() {
		List<List<ObservationReal>> newList = new ArrayList<>();
		if (train != null) {
			for (List<?> objects : train) {
				ArrayList<ObservationReal> tempList = (ArrayList<ObservationReal>) objects;
				newList.add(tempList);
			}
			if (k > 0 && numDist > 0) {
				kmeans = new KMeansLearner<>(k, new OpdfGaussianMixtureFactory(numDist), newList);
			}
		} else {
			//split the data into 1000bp chunks
			int a;
			int i;
			for (i = 0; i < obs.size() - 1000; i += 1000) {
				List<ObservationReal> temp = new ArrayList<>();
				for (a = i; a < i + 1000; a++) {
					ObservationReal o = (ObservationReal) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		}
		return newList;
	}
	
	/**
	 * Format the data as double data
	 *
	 * @return a List of Lists of ObservationReal representing the data for baum welch training
	 */
	@SuppressWarnings("unchecked")
	private List<List<ObservationReal>> formatRealTrain() {
		List<List<ObservationReal>> newList = new ArrayList<>();
		if (train != null) {
			for (List<?> objects : train) {
				ArrayList<ObservationReal> tempList = (ArrayList<ObservationReal>) objects;
				newList.add(tempList);
			}
			if (k > 0) {
				kmeans = new KMeansLearner<>(k, new OpdfGaussianFactory(), newList);
			}
		} else {
			//split the data into 1000bp chunks
			int a;
			int i;
			for (i = 0; i < obs.size() - 1000; i += 1000) {
				List<ObservationReal> temp = new ArrayList<>();
				for (a = i; a < i + 1000; a++) {
					ObservationReal o = (ObservationReal) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		}
		return newList;
	}
	
	/**
	 * Format the data as double data
	 *
	 * @return a List of ObservationReal representing the data for viterbi decoding
	 */
	@SuppressWarnings("unchecked")
	private List<ObservationReal> formatRealObs() {
		return (List<ObservationReal>) (ArrayList<ObservationReal>) obs;
	}
	
	/**
	 * Format hmm for double data
	 *
	 * @return a hmm
	 */
	@SuppressWarnings("unchecked")
	private Hmm<ObservationReal> formatRealHmm() {
		if (hmm != null) {
			return (Hmm<ObservationReal>) hmm;
		} else {
			return null;
		}
	}
	
	/**
	 * Format the data as integer data
	 *
	 * @return a List of List of ObservationInteger representing the data for baum welch training
	 */
	@SuppressWarnings("unchecked")
	private List<List<ObservationInteger>> formatIntTrain() {
		List<List<ObservationInteger>> newList = new ArrayList<>();
		if (train != null) {
			for (List<?> objects : train) {
				
				ArrayList<ObservationInteger> tempList = (ArrayList<ObservationInteger>) objects;
				newList.add(tempList);
			}
			if (k > 0 && numDist > 0) {
				kmeans = new KMeansLearner<>(k, new OpdfIntegerFactory(numDist), newList);
			}
		} else {
			//split the data into 1000bp chunks
			int a;
			int i;
			for (i = 0; i < obs.size() - 1000; i += 1000) {
				List<ObservationInteger> temp = new ArrayList<>();
				for (a = i; a < i + 1000; a++) {
					ObservationInteger o = (ObservationInteger) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		}
		return newList;
	}
	
	/**
	 * Format hmm for integer data
	 *
	 * @return a hmm
	 */
	@SuppressWarnings("unchecked")
	private Hmm<ObservationInteger> formatIntHmm() {
		if (hmm != null) {
			return (Hmm<ObservationInteger>) hmm;
		} else {
			return null;
		}
	}
	
	/**
	 * Format the data as integer data
	 *
	 * @return a List of ObservationInteger representing the data for viterbi decoding
	 */
	@SuppressWarnings("unchecked")
	private List<ObservationInteger> formatIntObs() {
		return (List<ObservationInteger>) (ArrayList<ObservationInteger>) obs;
	}
	
	/**
	 * Access whether the data is in integer form
	 *
	 * @return a boolean to determine if data is in integer form
	 */
	private boolean isInt() {
		return method.equals("Integer");
	}
	
	/**
	 * Access whether the data is in double form
	 *
	 * @return a boolean to determine if data is in double form
	 */
	private boolean isReal() {
		return method.equals("Real");
	}
	
	/**
	 * Access whether the data is in multivariate gaussian form
	 *
	 * @return a boolean to determine if data is in multivariate gaussian form
	 */
	private boolean isVector() {
		return method.equals("Vector");
	}
	
	/**
	 * Access whether the data is in mixture form
	 *
	 * @return a boolean to determine if data is in mixture form
	 */
	private boolean isMixture() {
		return method.equals("Mixture");
	}
	
	/**
	 * Access whether the model is valid
	 *
	 * @param h a hmm
	 * @return a boolean to determine if the model is valid
	 */
	private boolean bwConverge(Hmm<?> h) {
		double pi = h.getPi(0);
		return Double.isNaN(pi);
	}
}
