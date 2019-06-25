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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import net.sf.javaml.core.DenseInstance;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;

public class RandomInitHMM {
	//class variables
	private String input;
	private int numStates;
	private Hmm<?> hmm;
	
	//Main method variables
	private static int states = 0;
	private static String file = null;
	private static String output = null;
	
	public RandomInitHMM(String i,int n) throws FileNotFoundException{
		input = i;
		numStates = n;
		read();
	}
	public Hmm<?> getHMM(){return hmm;}

	private void read() throws FileNotFoundException{
		ObservationReader reader = new ObservationReader(input);
		String method = reader.getMethod();
		if (method.equals("Vector")){
			hmm = readVector(reader);
		}
		else if(method.equals("Real") ){
			hmm = readReal(reader);
		}
		
	}
	private void read2() throws FileNotFoundException{
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		ArrayList<double[]> data = new ArrayList<double[]>();
		int numFeat = 0;
		while (inFile.hasNext()){
			String line = inFile.nextLine();
			String[] features = line.split(",");
			double[] values = new double[features.length];
			numFeat = features.length;
			for (int i = 0;i < features.length;i++){
				values[i] = Double.parseDouble(features[i]);
			}
			data.add(values);
		}
		double[][] values = new double[data.size()][numFeat];
		for(int i = 0;i < data.size();i++){
			double[] temp = data.get(i);
			for (int a = 0;a < temp.length;a++){
				values[i][a] = temp[a];
			}
		}
		data = null;
		double[] mu = new double[numFeat];
		double[] var = new double[numFeat];
		Mean mean = new Mean();
		Variance variance = new Variance();
		for (int i = 0; i < numFeat;i++){
			mu[i] = mean.evaluate(values[i]);
			var[i] = variance.evaluate(values[i]);
		}
		List<OpdfMultiGaussian> opdf = new ArrayList<OpdfMultiGaussian>();
		for (int i = 0;i < numStates;i++){
			
		}
	}
	
	@SuppressWarnings({"unchecked"})
	private Hmm<ObservationVector> readVector(ObservationReader reader){
		
		ArrayList<ObservationVector> obs = (ArrayList<ObservationVector>) reader.getObs();
		double[] initial = setInitial();
		double[][] trans = setTrans();
		List<OpdfMultiGaussian> opdf = setVectorPDF(obs);
		
		Hmm<ObservationVector> h = new Hmm<ObservationVector>(initial, trans, opdf);
		
		
		return h;
		
	}
	@SuppressWarnings({ "unused", "unchecked" })
	private Hmm<ObservationReal> readReal(ObservationReader reader){
		ArrayList<ObservationReal> obs = (ArrayList<ObservationReal>) reader.getObs();
		double[] initial = setInitial();
		double[][] trans = setTrans();
		//TODO: write opdf maker for monovariate gaussian
		return null;
	}
	private List<OpdfGaussian> setRealPDF(ArrayList<ObservationReal> obs){
		//TODO: write this method
		return null;
	}
	
	private List<OpdfMultiGaussian> setVectorPDF(ArrayList<ObservationVector> obs){
		List<OpdfMultiGaussian> opdf = new ArrayList<OpdfMultiGaussian>();
		ObservationVector o = obs.get(0);
		int dim = o.dimension();
		double[] means = new double[dim];
		Mean mu = new Mean();
		Variance v = new Variance();
		double[] vars = new double[dim];
		double[] vals = new double[obs.size()];
		for (int i = 0;i < dim;i++){
			
			for (int a = 0;a < obs.size();a++){
				vals[a] = obs.get(a).value(i);
			}
			means[i] = mu.evaluate(vals, 0, vals.length);
			vars[i] = v.evaluate(vals);
		}
		double[][] Means = new double[numStates][dim];
		double[][] Vars = new double[numStates][dim];
		for (int i = 0; i < means.length;i++){
			double mean = means[i];
			double var = vars[i];
			
			double sd = Math.sqrt(var);
			double meanLower = mean - (3.0*sd);
			double meanUpper = mean + (3.0*sd);
			double varLower = 0.5 * var;
			double varUpper = 3 * var;
			
			
			double meanStep = (meanUpper - meanLower) / ((double)numStates - 1);
			double varStep = (varUpper - varLower) / ((double) numStates - 1);
			//System.out.println((meanUpper-meanLower)+"\t"+meanStep);
			for(int a = 0;a < numStates;a++){
				Means[a][i] = meanLower + (a * meanStep);
				Vars[a][i] = varLower + (a * varStep);
				//System.out.println(Means[a][i]);
				//System.out.println(Vars[a][i]);
			}
			
		}
		for (int i = 0;i < Means.length;i++){
			//System.out.println(Means.length);
			double[][] cov = new double[dim][dim];
			for (int a = 0;a < Means[i].length;a++){
				//System.out.println(Means[i].length);
				cov[a][a] = Vars[i][a];
				//System.out.println(cov[a][a]);
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(Means[i],cov);
			opdf.add(pdf);
		}
		
		return opdf;
		
	}
	private double[][] setTrans(){
		double[][] trans = new double[numStates][numStates];
		for (int i = 0;i < numStates;i++){
			for(int a = 0;a < numStates;a++){
				trans[i][a] = (double) 1/numStates;
			}
		}
		return trans;
	}
	private double[] setInitial(){
		double[] initial = new double[numStates];
		for(int i = 0;i < numStates;i++){
			initial[i] = (double) 1/numStates;
		}
		return initial;
	}
	
	public static void main(String[] args) throws IOException{
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			
			case'i':
				file = (args[i+1]);
				i++;
				break;
			case 'n':
				states = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'o':
				output = args[i+1];
				i++;
				break;
			
			}
		}
		if (file == null || states == 0 || output == null){
			printUsage();
			System.exit(1);
		}
		
		RandomInitHMM init = new RandomInitHMM(file,states);
		Hmm<?> hmm = init.getHMM();
		System.out.println(hmm.toString());
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, hmm);
		
	}
	private static void printUsage(){
		System.out.println("Usage: java -jar RandomInitHMM.jar");
		System.out.println("Required Parameters:");
		System.out.println("-i <File> Observation File in proper format");
		System.out.println("-n <int> Number of States");
		System.out.println("-o <File> Output file to write binary HMM");
	}
}
