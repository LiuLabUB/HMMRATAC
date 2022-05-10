package GEMM;

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

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class BasicEM {
	
	private final double[] mode;
	private final double[] weights;
	private final double[] mu;
	private final double[] lambda;
	private final double[] data;
	private int maxIter = 10;
	private double epsilon = 0.0005;
	private int numIter = 0;
	
	public BasicEM(double[] mode, double[] weights, double[] means, double[] lambda, double[] data) {
		this.mode = mode;
		this.weights = weights;
		this.mu = means;
		this.lambda = lambda;
		this.data = data;
	}
	
	public void setMaxIter(int m) {
		maxIter = m;
	}
	
	public void setEpsilon(double e) {
		epsilon = e;
	}
	
	public int getNumberIterations() {
		return numIter;
	}
	
	public double[] getMeans() {
		return mu;
	}
	
	public double[] getWeights() {
		return weights;
	}
	
	public double[] getLambda() {
		return lambda;
	}
	
	public void learn() {
		
		int i;
		double[] oldMu = new double[mode.length];
		double[] oldWeights = new double[mode.length];
		double[] oldLam = new double[mode.length];
		
		for (i = 0; i < maxIter; i++) {
			
			System.arraycopy(mu, 0, oldMu, 0, mode.length);
			System.arraycopy(lambda, 0, oldLam, 0, mode.length);
			System.arraycopy(weights, 0, oldWeights, 0, mode.length);
			
			iterate();
			
			int counter = 0;
			for (int a = 0; a < mode.length; a++) {
				if (converges(oldMu[a], mu[a], epsilon)
						&& converges(oldWeights[a], weights[a], epsilon)
						&& converges(oldLam[a], lambda[a], epsilon)) {
					counter++;
				}
			}
			if (counter == mode.length) {
				break;
			}
		}
		numIter = i;
	}
	
	private boolean converges(double value1, double value2, double epsilon) {
		return Math.abs(value1 - value2) <= epsilon;
	}
	
	public void iterate() {
		ArrayList<ArrayList<Double>> separatedData = new ArrayList<>();
		for (int i = 0; i < mode.length; i++)
			separatedData.add(new ArrayList<>());
		double[] temp = new double[mode.length];
		double counter = 0.0;
		for (double d : data) {
			for (int l = 0; l < mode.length; l++) {
				temp[l] = getWeightedDensity(mode[l], d, mu[l], lambda[l], weights[l]);
			}
			int index = returnGreater(temp);
			if (index != -1) {
				separatedData.get(index).add(d);
				counter++;
			}
		}
		
		for (int l = 0; l < mode.length; l++) {
			double[] d = convertFromList(separatedData.get(l));
			double jump = 1.5;
			if (mode[l] == 1) {
				mu[l] = getArithmeticMedian(d);
				lambda[l] = updateLaplaceBeta(d);
			} else if (mode[l] == 2) {
				mu[l] = mu[l] + (jump * (getArithmeticMean(d) - mu[l]));
				lambda[l] = lambda[l] + (jump * (getArithmeticStdDev(d) - lambda[l]));
			} else if (mode[l] == 0.5) {
				mu[l] = getArithmeticMean(d);
				lambda[l] = lambda[l];
			} else if (mode[l] == 0) {
				mu[l] = mu[l];
				lambda[l] = updateLaplaceBeta(d);
			} else if (mode[l] == 3) {
				mu[l] = mu[l];
				lambda[l] = lambda[l];
			}
			double xml = (double) d.length / counter;
			weights[l] = weights[l] + (jump * (xml - weights[l]));
		}
	}
	
	/*
	 * returns the index of data that is the greatest value or -1 if the largest value in data is duplicated
	 */
	private int returnGreater(double[] data) {
		int largest = -1;
		double greatest = -Double.MAX_VALUE;
		for (int i = 0; i < data.length; i++) {
			if (data[i] == greatest) {
				largest = -1;
			} else if (data[i] > greatest) {
				greatest = data[i];
				largest = i;
			}
		}
		return largest;
	}
	
	private double getWeightedDensity(double p, double x, double mean, double lambda, double weight) {
		return weight * getDensity(p, x, mean, lambda);
	}
	
	private double getDensity(double p, double x, double mean, double lambda) {
		double density = 0.0;
		
		if (p == 2) {
			NormalDistribution dis = new NormalDistribution(mean, lambda);
			density = dis.density(x);
		} else if (p == 0.5) {
			ExponentialDistribution dis = new ExponentialDistribution(mean);
			density = dis.density(x);
		}
		return density;
	}
	
	private static double[] convertFromList(ArrayList<Double> data) {
		return data.stream().mapToDouble(i -> i).toArray();
	}
	
	private double getArithmeticMean(double[] data) {
		return new Mean().evaluate(data);
	}
	
	private double getArithmeticStdDev(double[] data) {
		return new StandardDeviation().evaluate(data);
	}
	
	private double getArithmeticMedian(double[] data) {
		return new Percentile().evaluate(data, 50.0);
	}
	
	private double updateLaplaceBeta(double[] data) {
		double mean = getArithmeticMedian(data);
		return Arrays.stream(data).map(d -> Math.abs(d - mean)).sum() / data.length;
	}
	
	public static void main(String[] args) {
		Scanner inFile = null;
		try {
			inFile = new Scanner(new FileReader(args[0]));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<Double> temp = new ArrayList<>();
		while (inFile.hasNextLine()) {
			String line = inFile.nextLine();
			String[] feat = line.split("//s+");
			double value = Double.parseDouble(feat[0]);
			temp.add(value);
		}
		double[] data = convertFromList(temp);
		
		//idea is to use the local max and local mins to identify the initial parameters. short initial mean will be readlength
		//mono initial could be first local max (certain distance from short) and di/tri could be 2x or 3x of mono etc
		
		double[] weights = new double[3];
		double[] means = new double[3];
		double[] lam = new double[3];
		double[] mode = new double[3];
		
		means[0] = 150;
		means[1] = 300;
		means[2] = 450;
		lam[0] = lam[1] = lam[2] = 20;
		mode[0] = mode[1] = mode[2] = 2;
		
		double cut1 = (means[1] - means[0]) / 2 + means[0];
		double cut2 = (means[2] - means[1]) / 2 + means[1];
		int counter1 = 0;
		int counter2 = 0;
		int counter3 = 0;
		for (double d : data) {
			if (d < cut1)
				counter1++;
			else if (d >= cut1 && d <= cut2)
				counter2++;
			else
				counter3++;
		}
		weights[0] = (double) counter1 / data.length;
		weights[1] = (double) counter2 / data.length;
		weights[2] = (double) counter3 / data.length;
		
		BasicEM em = new BasicEM(mode, weights, means, lam, data);
		em.setMaxIter(50);
		em.learn();
	}
	
	public static int[] getInflectionPoints(double[] smoothed, int order) {
		int[] local;
		ArrayList<Integer> temp = new ArrayList<>();
		for (int i = order; i < smoothed.length - order; i++) {
			double upSlope = Math.abs(smoothed[i - order] - smoothed[i]) / (double) order + 1.0;
			double downSlope = Math.abs(smoothed[i] - smoothed[i + order]) / (double) order + 1.0;
			if (upSlope != downSlope) {
				temp.add(i);
			}
		}
		local = new int[temp.size()];
		for (int i = 0; i < temp.size(); i++)
			local[i] = temp.get(i);
		return local;
	}
	
	public static int[] getLocalExtrema(double[] smoothed, String comparator, int order) {
		int[] local;
		ArrayList<Integer> temp = new ArrayList<>();
		for (int i = (order); i < smoothed.length - (order); i++) {
			if (comparator.toLowerCase().contains("great")) {
				if (Math.max(smoothed[i], smoothed[i - order]) == smoothed[i] &&
						Math.max(smoothed[i], smoothed[i + order]) == smoothed[i] &&
						smoothed[i] != smoothed[i + order] &&
						smoothed[i] != smoothed[i - order]) {
					temp.add(i);
				}
			} else if (comparator.toLowerCase().contains("less")) {
				if (Math.min(smoothed[i], smoothed[i - order]) == smoothed[i] &&
						Math.min(smoothed[i], smoothed[i + order]) == smoothed[i] &&
						smoothed[i] != smoothed[i + order] &&
						smoothed[i] != smoothed[i - order]) {
					temp.add(i);
				}
			}
		}
		local = new int[temp.size()];
		for (int i = 0; i < temp.size(); i++) {
			local[i] = temp.get(i);
		}
		return local;
	}
	
	public static int[] getLocalMaxMin(double[] _smoothed) {
		int[] local = new int[_smoothed.length];
		for (int i = 2; i < _smoothed.length - 2; i++) {
			//Calculates all local maxima and local minima points in the smoothed pileups
			
			double pointone = _smoothed[i - 2];
			double pointtwo = _smoothed[i - 1];
			double pointthree = _smoothed[i];
			double pointfour = _smoothed[i + 1];
			double pointfive = _smoothed[i + 2];
			
			double slopepointtwo = (pointthree - pointone) / ((i) - (i - 2));
			double slopepointthree = (pointfour - pointtwo) / ((i + 1) - (i - 1));
			double slopepointfour = (pointfive - pointthree) / ((i + 2) - (i));
			
			if (slopepointthree == 0) {
				double secondslopepointthree = (slopepointfour - slopepointtwo) / ((i + 1) - (i - 1));
				if (secondslopepointthree > 0) {
					local[i]--;
				}
				if (secondslopepointthree < 0) {
					local[i]++;
				}
			}
		}
		return local;
	}
}
