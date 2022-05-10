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

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class SlidingViterbiDriver {
	
	private static String obs = null;
	private static String train = null;
	private static File hmm = null;
	private static boolean welch = false;
	private static int k = 0;
	private static int numDist = 0;
	private static int winSize = 1000;
	private static int slideSize = 1000;
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {
			
			switch (args[i].charAt((1))) {
				
				case 'o':
					obs = (args[i + 1]);
					i++;
					break;
				case 't':
					train = (args[i + 1]);
					i++;
					break;
				case 'm':
					hmm = new File(args[i + 1]);
					i++;
					break;
				case 'w':
					String w = args[i + 1].toLowerCase();
					if (w.equals("t") || w.equals("true")) {
						welch = true;
					}
					i++;
					break;
				case 'k':
					k = Integer.parseInt(args[i + 1]);
					i++;
					break;
				case 'n':
					numDist = Integer.parseInt(args[i + 1]);
					i++;
					break;
				case 's':
					winSize = Integer.parseInt(args[i + 1]);
					i++;
					break;
				case 'l':
					slideSize = Integer.parseInt(args[i + 1]);
					i++;
					break;
				case 'h':
					printUsage();
					System.exit(0);
			}
		}
		if (obs == null || (train == null && hmm == null)
				|| (train == null && welch)
				|| (train == null && k > 0)) {
			printUsage();
			System.exit(0);
		}
		
		
		TrainingReader trainReader;
		List<List<?>> trainList = null;
		if (train != null) {
			trainReader = new TrainingReader(train);
			trainList = trainReader.getObs();
		}
		Hmm<?> h = null;
		if (hmm != null) {
			h = HmmBinaryReader.read(new FileInputStream(hmm));
		}
		ObservationReader obsReader = new ObservationReader(obs);
		List<?> obsList = obsReader.getObs();
		String method = obsReader.getMethod();
		
		ArrayList<Points> results = new ArrayList<>();
		int i;
		for (i = 0; i < obsList.size() - winSize; i += slideSize) {
			
			List<?> temp = run(obsList, method, winSize, i);
			runViterbi(trainList, h, method, results, i, temp, winSize);
			
		}
		winSize = obsList.size() - i;
		List<?> temp = run(obsList, method, winSize, i);
		runViterbi(trainList, h, method, results, i, temp, winSize);
		
		
		double[] initial = new double[results.size()];
		for (int a = 0; a < results.size(); a++) {
			initial[a] = results.get(a).getProb();
		}
		
		for (Points result : results) {
			List<?> tempList = run(obsList, method, result.getStop(),
					result.getStart());
			ResultNode node = runViterbi2(trainList, h, method, tempList);
			int[] states = node.getStates();
			for (int x = 0; x < states.length; x++) {
				String obs = tempList.get(x).toString();
				System.out.println(result.getStart() + x + "\t" + states[x] + "\t" + result.getProb() + "\t" + obs);
			}
		}
	}
	
	
	private static List<?> run(List<?> obsList, String method, int winSize, int i) {
		if (Objects.equals(method, "Int")) {
			List<ObservationInteger> tempList = new ArrayList<>();
			for (int a = i; a < i + winSize; a++) {
				ObservationInteger o = (ObservationInteger) obsList.get(a);
				tempList.add(o);
			}
			return tempList;
		} else if (Objects.equals(method, "Real") || Objects.equals(method, "Mixture")) {
			List<ObservationReal> tempList = new ArrayList<>();
			for (int a = i; a < i + winSize; a++) {
				ObservationReal o = (ObservationReal) obsList.get(a);
				tempList.add(o);
			}
			return tempList;
		} else if (Objects.equals(method, "Vector")) {
			List<ObservationVector> tempList = new ArrayList<>();
			for (int a = i; a < i + winSize; a++) {
				ObservationVector o = (ObservationVector) obsList.get(a);
				tempList.add(o);
			}
			return tempList;
		} else {
			return null;
		}
	}
	
	
	private static void runViterbi(List<List<?>> trainList, Hmm<?> h,
								   String method, ArrayList<Points> results, int i,
								   List<?> tempList, int winSize) {
		RobustHMM HMM = new RobustHMM(tempList, trainList, h, welch, k, method, numDist);
		double prob = HMM.getProb2();
		Points p = new Points(prob, i, winSize);
		results.add(p);
	}
	
	private static ResultNode runViterbi2(List<List<?>> trainList, Hmm<?> h, String method, List<?> tempList) {
		RobustHMM HMM = new RobustHMM(tempList, trainList, h, welch, k, method, numDist);
		double prob = HMM.getProb();
		int[] states = HMM.getStates();
		return new ResultNode(prob, states);
	}
	
	
	private static void printUsage() {
		System.out.println("Usage: java -jar ViterbiDriver.jar");
		System.out.println("Required Parameters:");
		System.out.println("-o <File> Observation File in specific format");
		System.out.println("Optional Parameters:");
		System.out.println("-t <File> Training set file in specific format. Note if -w is true or -h is missing or -k is set, this is required");
		System.out.println("-m <File> File containing binary representation of HMM. If missing, -t is required");
		System.out.println("-w <T/F> Whether or not to use Baum Welch. If true, -t is required.");
		System.out.println("-k <int> Number of states used for KMeans training. If set, -t is required and -h is not required");
		System.out.println("-n <int> Number of Gaussian Distributions if using a mixture model or the number of discrete integers used for Integer model. Note: required if -k is set");
		System.out.println("-s <int> Size of sliding window. Default = 1000");
		System.out.println("-l <int> Size of the slide. Default = 1000");
		System.out.println("-h Print this help message and exit.");
	}
}

class Points {
	private double prob;
	private int start;
	private int stop;
	
	public Points(double pr, int s, int st) {
		prob = pr;
		start = s;
		stop = st;
	}
	
	public double getProb() {
		return prob;
	}
	
	public void setProb(double p) {
		prob = p;
	}
	
	public int getStart() {
		return start;
	}
	
	public void setStart(int s) {
		start = s;
	}
	
	public int getStop() {
		return stop;
	}
	
	public void setStop(int s) {
		stop = s;
	}
}

class ResultNode {
	private double prob;
	private int[] states;
	
	public ResultNode(double p, int[] s) {
		prob = p;
		states = s;
	}
	
	public double getProb() {
		return prob;
	}
	
	public void setProb(double p) {
		prob = p;
	}
	
	public int[] getStates() {
		return states;
	}
	
	public void setStates(int[] s) {
		states = s;
	}
}
