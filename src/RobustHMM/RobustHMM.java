package RobustHMM;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixtureFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfIntegerFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;

public class RobustHMM {
	//inputs
	private List<List<?>> train;
	private List<?> obs;
	private Hmm<?> hmm;
	private boolean welch;
	private int k;
	private String method;
	private int numDist;
	private int BWIter = 10;
	//specifics
	private KMeansLearner<?> kmeans;
	private BaumWelchLearner bw;
	private BaumWelchScaledLearner sbw;
	
	private int[] states;
	private Hmm<?> outHmm;
	private double prob;
	private double prob2;
	
	public RobustHMM(List<?> o,List<List<?>> t,Hmm<?> h,boolean w,int K,String m,int n){
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
	public void setBaumWelchMaxIter(int n){
		BWIter = n;
	}
	public int[] getStates(){
		return states;
	}
	public double getProb(){
		return prob;
	}
	public double getProb2(){return prob2;}
	
	public Hmm<?> getHmm(){return outHmm;}
	private void run(){
		if (isInt()){
			runInt();
		}
		else if (isReal()){
			runReal();
		}
		else if (isVector()){
			runVector();
		}
		else if (isMixture()){
			runMixture();
		}
	}
	@SuppressWarnings("unchecked")
	private void runMixture(){
		List<List<ObservationReal>> newTrain = formatMixtureTrain();
		List<ObservationReal> newObs = formatRealObs();
		Hmm<ObservationReal> newHmm = null;
		if (k > 0){
			newHmm = (Hmm<ObservationReal>) kmeans.learn();
		}
		else{
			newHmm = formatRealHmm();
		}
		if (welch == false){
			states = viterbiReal(newHmm,newObs);
			outHmm = newHmm;
		}
		else{
			for (int x = 0;x < BWIter;x++){
				newHmm = bw.iterate(newHmm, newTrain);
				if(bwConverge(newHmm)){
					break;
				}
			}
			if (!bwConverge(newHmm)){
				for (int y = 0;y < BWIter;y++){
					newHmm = sbw.iterate(newHmm, newTrain);
					if(bwConverge(newHmm)){
						break;
					}
				}
			}
			if (!bwConverge(newHmm)){
				System.out.println("Baum-Welch Iterations"+"\t"+BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
			states = viterbiReal(newHmm,newObs);
			outHmm = newHmm;
		}
	}
	@SuppressWarnings("unchecked")
	private void runVector(){
		List<List<ObservationVector>> newTrain = formatVectorTrain();
		List<ObservationVector> newObs = formatVectorObs();
		Hmm<ObservationVector> newHmm = null;
		if (k > 0){
			newHmm = (Hmm<ObservationVector>) kmeans.learn();
		}
		else{
			newHmm = formatVectorHmm();
		}
		if (welch == false){
			states = viterbiVector(newHmm,newObs);
			outHmm = newHmm;
		}
		else{
			
			for (int x = 0;x < BWIter;x++){
				newHmm = sbw.iterate(newHmm, newTrain);
				if(bwConverge(newHmm)){
					System.out.println("Baum Welch Error!!!");
					System.exit(1);
					break;
				}
			}
				
			/*
			if (!bwConverge(newHmm)){
				for (int y = 0;y < BWIter;y++){
					newHmm = sbw.iterate(newHmm, newTrain);
					if(bwConverge(newHmm)){
						break;
					}
				}
			}
			if (!bwConverge(newHmm)){
				System.out.println("Baum-Welch Iterations"+"\t"+BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}*/
			states = viterbiVector(newHmm,newObs);
			outHmm = newHmm;
		}
	}
	@SuppressWarnings("unchecked")
	private void runReal(){
		List<List<ObservationReal>> newTrain = formatRealTrain();
		List<ObservationReal> newObs = formatRealObs();
		Hmm<ObservationReal> newHmm = null;
		if (k > 0){
			newHmm = (Hmm<ObservationReal>) kmeans.learn();
		}
		else{
			newHmm = formatRealHmm();
		}
		if (welch == false){
			states = viterbiReal(newHmm,newObs);
			outHmm = newHmm;
		}
		else{
			for (int x = 0;x < BWIter;x++){
				newHmm = bw.iterate(newHmm, newTrain);
				if(bwConverge(newHmm)){
					break;
				}
			}
			if (!bwConverge(newHmm)){
				for (int y = 0;y < BWIter;y++){
					newHmm = sbw.iterate(newHmm, newTrain);
					if(bwConverge(newHmm)){
						break;
					}
				}
			}
			if (!bwConverge(newHmm)){
				System.out.println("Baum-Welch Iterations"+"\t"+BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
			states = viterbiReal(newHmm,newObs);
			outHmm = newHmm;
		}
	}
	@SuppressWarnings("unchecked")
	private void runInt(){
		List<List<ObservationInteger>> newTrain = formatIntTrain();
		List<ObservationInteger> newObs = formatIntObs();
		Hmm<ObservationInteger> newHmm = null;
		if (k > 0){
			newHmm = (Hmm<ObservationInteger>) kmeans.learn();
		}
		else{
			newHmm = formatIntHmm();
		}
		if (welch == false){
			states = viterbiInt(newHmm,newObs);
			outHmm = newHmm;
		}
		else{
			for (int x = 0;x < BWIter;x++){
				newHmm = bw.iterate(newHmm, newTrain);
				if(bwConverge(newHmm)){
					break;
				}
			}
			if (!bwConverge(newHmm)){
				for (int y = 0;y < BWIter;y++){
					newHmm = sbw.iterate(newHmm, newTrain);
					if(bwConverge(newHmm)){
						break;
					}
				}
			}
			if (!bwConverge(newHmm)){
				System.out.println("Baum-Welch Iterations"+"\t"+BWIter);
				throw new BaumWelchException("Baum Welch Failed!! Check your training set");
			}
			states = viterbiInt(newHmm,newObs);
			outHmm = newHmm;
		}
	}
	private int[] viterbiInt(Hmm<ObservationInteger> h,List<ObservationInteger> l){
		ViterbiCalculator vit = new ViterbiCalculator(l,h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	private int[] viterbiReal(Hmm<ObservationReal> h,List<ObservationReal> l){
		ViterbiCalculator vit = new ViterbiCalculator(l,h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	private int[] viterbiVector(Hmm<ObservationVector> h,List<ObservationVector> l){
		ViterbiCalculator vit = new ViterbiCalculator(l,h);
		int[] states = vit.stateSequence();
		prob = vit.lnProbability();
		prob2 = h.probability(l);
		return states;
	}
	
	@SuppressWarnings("unchecked")
	private List<List<ObservationVector>> formatVectorTrain(){
		if (train != null){
			List<List<ObservationVector>> newList = new ArrayList<List<ObservationVector>>();
			for (int i = 0;i < train.size();i++){
				ArrayList<ObservationVector> tempList = (ArrayList<ObservationVector>) train.get(i);
				newList.add(tempList);
			}
			int dim = newList.get(0).get(0).dimension();
			if (k > 0){
				kmeans = new KMeansLearner<ObservationVector>(k,new OpdfMultiGaussianFactory(dim),newList);
			}
			return newList;
		}
		/*
		 * This is a test. If no training set is given, duplicate the observation set and use that. 
		 * Could also attempt to split up the observation set into smaller, random pieces.
		 * 	Duplicating observation set resulted in out of memory error
		 * 	Spliting the data into random chunks resulted in non convergence of BW
		 * Try splitting the data into 1000bp chunks (for whole data)
		 * Note: either of these approaches may not work. 
		 * Update 8/27/15:
		 * 	Splitting into 1000bp chunks sometimes worked. Problem was that although individual tracks were not
		 * 	correlated in whole dataset, they sometimes were within the chunks, leading to problems
		 * 	Try splitting whole dataset in half instead. Less likely for data tracks to be correlated
		 * Original approach was simply to return null if train == null
		 */
		else{
			List<List<ObservationVector>> newList = new ArrayList<List<ObservationVector>>();
			//split the data into 1000bp chunks
			int a; int i;
			int halfSize = obs.size()/2;
			//System.out.println(halfSize);
			for (i = 0;i < obs.size()-1;i+=halfSize){
				//System.out.println("i is "+"\t"+i);
				List<ObservationVector> temp = new ArrayList<ObservationVector>();
				for (a = i;a < i+halfSize-1;a++){
					//System.out.println("A is "+"\t"+a);
					ObservationVector o = (ObservationVector) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
			
			
			return newList;
		}
	}
	private boolean inBetween(int start,int stop,int value){
		if ((value >= start && value <= stop) ||  ((value+15000) >= start && (value+15000)<=stop)){
			return true;
		}
		else{return false;}
	}
	
	@SuppressWarnings("unchecked")
	private List<ObservationVector> formatVectorObs(){
		List<ObservationVector> newObs = (ArrayList<ObservationVector>) obs;
		return newObs;
	}
	@SuppressWarnings("unchecked")
	private Hmm<ObservationVector> formatVectorHmm(){
		if (hmm != null){
		Hmm<ObservationVector> newHmm = (Hmm<ObservationVector>) hmm;
		return newHmm;}
		else {return null;}
	}
	@SuppressWarnings("unchecked")
	private List<List<ObservationReal>> formatMixtureTrain(){
		if (train != null){
			List<List<ObservationReal>> newList = new ArrayList<List<ObservationReal>>();
			for (int i = 0;i < train.size();i++){
				ArrayList<ObservationReal> tempList = (ArrayList<ObservationReal>) train.get(i);
				newList.add(tempList);
			}
			if (k >0 && numDist > 0){
				kmeans = new KMeansLearner<ObservationReal>(k,new OpdfGaussianMixtureFactory(numDist),newList);
			}
			return newList;}
		else{
			List<List<ObservationReal>> newList = new ArrayList<List<ObservationReal>>();
			//split the data into 1000bp chunks
			int a; int i;
			for (i = 0;i < obs.size()-1000;i+=1000){
				List<ObservationReal> temp = new ArrayList<ObservationReal>();
				for (a = i;a < i+1000;a++){
					ObservationReal o = (ObservationReal) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		
		
			return newList;
		}
	}
	
	@SuppressWarnings("unchecked")
	private List<List<ObservationReal>> formatRealTrain(){
		if (train != null){
		List<List<ObservationReal>> newList = new ArrayList<List<ObservationReal>>();
		for (int i = 0;i < train.size();i++){
			ArrayList<ObservationReal> tempList = (ArrayList<ObservationReal>) train.get(i);
			newList.add(tempList);
		}
		if (k > 0){
			kmeans = new KMeansLearner<ObservationReal>(k,new OpdfGaussianFactory(),newList);
		}
		return newList;}
		else{
			List<List<ObservationReal>> newList = new ArrayList<List<ObservationReal>>();
			//split the data into 1000bp chunks
			int a; int i;
			for (i = 0;i < obs.size()-1000;i+=1000){
				List<ObservationReal> temp = new ArrayList<ObservationReal>();
				for (a = i;a < i+1000;a++){
					ObservationReal o = (ObservationReal) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		
		
			return newList;
		}
	}
	@SuppressWarnings("unchecked")
	private List<ObservationReal> formatRealObs(){
		List<ObservationReal> newObs = (ArrayList<ObservationReal>) obs;
		return newObs;
	}
	@SuppressWarnings("unchecked")
	private Hmm<ObservationReal> formatRealHmm(){
		if (hmm!=null){
		Hmm<ObservationReal> newHmm = (Hmm<ObservationReal>) hmm;
		return newHmm;}
		else{return null;}
	}
	@SuppressWarnings("unchecked")
	private List<List<ObservationInteger>> formatIntTrain(){
		if (train != null){
			List<List<ObservationInteger>> newList = new ArrayList<List<ObservationInteger>>();
			for (int i = 0;i < train.size();i++){
			
				ArrayList<ObservationInteger> tempList = (ArrayList<ObservationInteger>) train.get(i);
				newList.add(tempList);
			}
			if (k > 0 && numDist > 0){
				kmeans = new KMeansLearner<ObservationInteger>(k,new OpdfIntegerFactory(numDist),newList);
			}
			return newList;}
		else{
			List<List<ObservationInteger>> newList = new ArrayList<List<ObservationInteger>>();
			//split the data into 1000bp chunks
			int a; int i;
			for (i = 0;i < obs.size()-1000;i+=1000){
				List<ObservationInteger> temp = new ArrayList<ObservationInteger>();
				for (a = i;a < i+1000;a++){
					ObservationInteger o = (ObservationInteger) obs.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		
		
			return newList;
		}
	}
	@SuppressWarnings("unchecked")
	private Hmm<ObservationInteger> formatIntHmm(){
		if (hmm!=null){
		Hmm<ObservationInteger> newHmm = (Hmm<ObservationInteger>) hmm;
		return newHmm;}
		else{return null;}
	}
	@SuppressWarnings("unchecked")
	private List<ObservationInteger> formatIntObs(){
		List<ObservationInteger> newObs = (ArrayList<ObservationInteger>) obs;
		return newObs;
	}
	private boolean isInt(){
		if (method.equals("Integer")){return true;}
		else{return false;}
	}
	private boolean isReal(){
		if (method.equals("Real")){return true;}
		else{return false;}
	}
	private boolean isVector(){
		if(method.equals("Vector")){return true;}
		else{return false;}
	}
	private boolean isMixture(){
		if(method.equals("Mixture")){return true;}
		else{return false;}
	}
	private boolean bwConverge(Hmm<?> h){
		double pi = h.getPi(0);
		if (Double.isNaN(pi)){
			return true;
		}
		else{
			return false;
		}
	}
}
