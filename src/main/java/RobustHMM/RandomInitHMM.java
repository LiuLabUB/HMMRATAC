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
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class RandomInitHMM {
	//class variables
	private final String input;
	private final int numStates;
	private Hmm<?> hmm;
	
	//Main method variables
	private static int states = 0;
	private static String file = null;
	private static String output = null;
	
	public RandomInitHMM(String i, int n) throws FileNotFoundException {
		input = i;
		numStates = n;
		read();
	}
	
	public Hmm<?> getHMM() {
		return hmm;
	}
	
	private void read() throws FileNotFoundException {
		ObservationReader reader = new ObservationReader(input);
		String method = reader.getMethod();
		if (method.equals("Vector")) {
			hmm = readVector(reader);
		} else if (method.equals("Real")) {
			hmm = readReal(reader);
		}
	}
	
	@SuppressWarnings({"unchecked"})
	private Hmm<ObservationVector> readVector(ObservationReader reader) {
		
		ArrayList<ObservationVector> obs = (ArrayList<ObservationVector>) reader.getObs();
		double[] initial = setInitial();
		double[][] trans = setTrans();
		List<OpdfMultiGaussian> opdf = setVectorPDF(obs);
		
		return new Hmm<>(initial, trans, opdf);
		
	}
	
	@SuppressWarnings({"unused", "unchecked"})
	private Hmm<ObservationReal> readReal(ObservationReader reader) {
		ArrayList<ObservationReal> obs = (ArrayList<ObservationReal>) reader.getObs();
		double[] initial = setInitial();
		double[][] trans = setTrans();
		//TODO: write opdf maker for monovariate gaussian
		return null;
	}
	
	private List<OpdfGaussian> setRealPDF(ArrayList<ObservationReal> obs) {
		//TODO: write this method
		return null;
	}
	
	private List<OpdfMultiGaussian> setVectorPDF(ArrayList<ObservationVector> obs) {
		List<OpdfMultiGaussian> opdf = new ArrayList<>();
		ObservationVector o = obs.get(0);
		int dim = o.dimension();
		double[] means = new double[dim];
		Mean mu = new Mean();
		Variance v = new Variance();
		double[] vars = new double[dim];
		double[] vals = new double[obs.size()];
		for (int i = 0; i < dim; i++) {
			
			for (int a = 0; a < obs.size(); a++) {
				vals[a] = obs.get(a).value(i);
			}
			means[i] = mu.evaluate(vals, 0, vals.length);
			vars[i] = v.evaluate(vals);
		}
		double[][] Means = new double[numStates][dim];
		double[][] Vars = new double[numStates][dim];
		for (int i = 0; i < means.length; i++) {
			double mean = means[i];
			double var = vars[i];
			
			double sd = Math.sqrt(var);
			double meanLower = mean - (3.0 * sd);
			double meanUpper = mean + (3.0 * sd);
			double varLower = 0.5 * var;
			double varUpper = 3 * var;
			
			
			double meanStep = (meanUpper - meanLower) / ((double) numStates - 1);
			double varStep = (varUpper - varLower) / ((double) numStates - 1);
			for (int a = 0; a < numStates; a++) {
				Means[a][i] = meanLower + (a * meanStep);
				Vars[a][i] = varLower + (a * varStep);
			}
		}
		for (int i = 0; i < Means.length; i++) {
			double[][] cov = new double[dim][dim];
			for (int a = 0; a < Means[i].length; a++) {
				cov[a][a] = Vars[i][a];
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(Means[i], cov);
			opdf.add(pdf);
		}
		
		return opdf;
		
	}
	
	private double[][] setTrans() {
		double[][] trans = new double[numStates][numStates];
		for (int i = 0; i < numStates; i++) {
			for (int a = 0; a < numStates; a++) {
				trans[i][a] = (double) 1 / numStates;
			}
		}
		return trans;
	}
	
	private double[] setInitial() {
		double[] initial = new double[numStates];
		for (int i = 0; i < numStates; i++) {
			initial[i] = (double) 1 / numStates;
		}
		return initial;
	}
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {
			
			switch (args[i].charAt((1))) {
				
				case 'i':
					file = (args[i + 1]);
					i++;
					break;
				case 'n':
					states = Integer.parseInt(args[i + 1]);
					i++;
					break;
				case 'o':
					output = args[i + 1];
					i++;
					break;
			}
		}
		if (file == null || states == 0 || output == null) {
			printUsage();
			System.exit(1);
		}
		
		RandomInitHMM init = new RandomInitHMM(file, states);
		Hmm<?> hmm = init.getHMM();
		System.out.println(hmm.toString());
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter.write(out, hmm);
	}
	
	private static void printUsage() {
		System.out.println("Usage: java -jar RandomInitHMM.jar");
		System.out.println("Required Parameters:");
		System.out.println("-i <File> Observation File in proper format");
		System.out.println("-n <int> Number of States");
		System.out.println("-o <File> Output file to write binary HMM");
	}
}
