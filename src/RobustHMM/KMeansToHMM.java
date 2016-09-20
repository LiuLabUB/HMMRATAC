package RobustHMM;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.math3.stat.correlation.StorelessCovariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import stats.KMeans2;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;

public class KMeansToHMM {

	private static String input = null;
	private static Hmm<?> hmm = null;
	private static int K = 4;
	private static String output = null;
	private static int numIter = 100;
	private static boolean diagCov = true;
	private static boolean equalInitial = true;
	private static boolean equalTrans = false;
	
	public static void main(String[] args) throws IOException{
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {

			case'i':
				input = args[i+1];
				i++;
				break;
			case'k':
				K = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'n':
				numIter = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'o':
				output = args[i+1];
				i++;
				break;
			case'd':
				String val = args[i+1].toLowerCase();
				if (val.contains("t")){
					diagCov = true;
				}
				else{diagCov = false;}
				i++;
				break;
			case'e':
				String value = args[i+1].toLowerCase();
				if (value.contains("t")){
					equalInitial = true;
				}
				else{equalInitial = false;}
				i++;
				break;
			case't':
				
				equalTrans = true;
				
				i++;
				break;
			case'h':
				printUsage();
				
				System.exit(1);
			}
		}
		
		if (input == null || output == null){
			printUsage();
			
			System.exit(1);
		}
		long start = System.nanoTime();
		KMeansToHMM convert = new KMeansToHMM();
		Hmm<?> hmm = convert.getHMM();
		System.out.println(hmm.toString());
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, hmm);
		long end = System.nanoTime();
		long length = end - start;
		long seconds = length / 1000000000;
		System.out.println("****Time For KMeans to HMM =\t"+seconds+" seconds****");
		
	}//main
	
	public static Hmm<?> averageModels(ArrayList<Hmm<?>> models){
		return null;
		//TODO: make method to average a list of multiple, ordered models
	}
	
	public KMeansToHMM() throws FileNotFoundException{
		build(input,K,numIter, diagCov, equalInitial,equalTrans);
		sort((Hmm<ObservationVector>) hmm);
	}
	public KMeansToHMM(Dataset d,int K,int numIter,boolean diag,boolean equal,boolean equal2){
		build2(d,K,numIter, diag, equal,equal2);
		sort((Hmm<ObservationVector>) hmm);
		
	}
	
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
	
	public Hmm<?> getHMM(){return hmm;}
	private void build2(Dataset data,int K,int numIter,boolean diag,boolean equal,boolean equal2){
		int numFeatures = data.noAttributes();
		//System.out.println("Number of features\t"+numFeatures);
		//KMeans2 kmeans = new KMeans2(data,K,numIter);
		KMeans k = new KMeans(K,numIter);
		Dataset[] clustered = k.cluster(data);
		//System.out.println("Kmeans done");
		int[] assignments = new int[data.size()];
		List<OpdfMultiGaussian> opdf = new ArrayList<OpdfMultiGaussian>();
		
		for (int a = 0; a < clustered.length;a++){
			//System.out.println(a);
			Dataset cluster = clustered[a];
			//System.out.println(cluster.size());
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
			//printCovariance(covMat);
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
	
	private void build(String input,int K,int numIter,boolean diag,boolean equal,boolean equal2) throws FileNotFoundException{
		int numFeatures = 0;
		Dataset data = new DefaultDataset();
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		while (inFile.hasNext()){
			String line = inFile.nextLine();
			String[] features = line.split(",");
			double[] values = new double[features.length];
			numFeatures = features.length;
			for (int i = 0;i < features.length;i++){
				values[i] = Double.parseDouble(features[i]);
			}
			DenseInstance ins = new DenseInstance(values);
			data.add(ins);
		}
		
		KMeans2 kmeans = new KMeans2(data,K,numIter);
		Dataset[] clustered = kmeans.cluster();
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
			printCovariance(covMat);
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
	public void printCovariance(double[][] cov){
		for (int i = 0;i < cov.length;i++){
			for (int a = 0;a < cov[i].length;a++){
				if(a < cov[i].length-1){
					System.out.print(cov[i][a]+"\t");
				}
				else{
					System.out.println(cov[i][a]);
				}
			}
		}
	}
	
	public static void printUsage(){
		System.out.println("Usage: java -jar KMeansToHMM.jar");
		System.out.println("Required Paramters:");
		System.out.println("-i <File> CSV file containing values to cluster");
		System.out.println("-o <File> Output file to write HMM to");
		System.out.println("Optional Paramters:");
		System.out.println("-k <int> Number of clusters. Default = 4");
		System.out.println("-n <int> Number of iterations of KMeans algorithm. Default = 100");
		System.out.println("-d <Boolean> Whether or not to make covariance matrices diagonal");
		System.out.println("-t <boolean> Whether or not to make the transition matrix proportional");
	}
	
}
