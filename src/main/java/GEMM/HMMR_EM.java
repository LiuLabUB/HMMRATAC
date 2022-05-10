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

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class HMMR_EM {
	
	private final double[] weights;
	private final double[] mu;
	private final double[] lambda;
	private final double[] data;
	
	/**
	 * Constructor for generating new HMMR_EM object
	 *
	 * @param weights an Array of doubles representing the weights of the distribution
	 * @param means an Array of doubles representing the means of the distribution
	 * @param lambda an Array of doubles representing the standard deviations of the distribution
	 * @param data an Array of doubles representing the read length data
	 */
	public HMMR_EM(double[] weights, double[] means, double[] lambda, double[] data) {
		this.weights = weights;
		this.mu = means;
		this.lambda = lambda;
		this.data = data;
	}
	
	/**
	 * Access the means
	 *
	 * @return an Array of doubles representing the means
	 */
	public double[] getMeans() {
		return mu;
	}
	
	/**
	 * Access the standard deviations
	 *
	 * @return an Array of doubles representing the standard deviations
	 */
	public double[] getLambda() {
		return lambda;
	}
	
	/**
	 * Perform one iteration of the EM algorithm
	 */
	public void iterate() {
		
		double[] temp = new double[mu.length];
		double[] counter = new double[mu.length];
		double total = 0.0;
		Mean[] means = new Mean[mu.length];
		StandardDeviation[] std = new StandardDeviation[mu.length];
		for (int i = 0; i < means.length; i++) {
			means[i] = new Mean();
			std[i] = new StandardDeviation();
		}
		for (double d : data) {
			for (int l = 0; l < mu.length; l++) {
				temp[l] = getWeightedDensity(d, mu[l], lambda[l], weights[l]);
			}
			int index = returnGreater(temp);
			if (index != -1) {
				means[index].increment(d);
				std[index].increment(d);
				counter[index]++;
				total++;
			}
		}
		for (int i = 0; i < mu.length; i++) {
			double jump = 1.5;
			mu[i] = mu[i] + (jump * (means[i].getResult() - mu[i]));
			lambda[i] = lambda[i] + (jump * (std[i].getResult() - lambda[i]));
			weights[i] = counter[i] / total;
		}
	}
	
	/**
	 * Iterate through the EM algorithm until data convergence
	 */
	public void learn() {
		
		double[] oldMu = new double[mu.length];
		double[] oldWeights = new double[mu.length];
		double[] oldLam = new double[mu.length];
		boolean converged = false;
		int iter = 0;
		while (!converged) {
			
			for (int a = 0; a < mu.length; a++) {
				oldMu[a] = mu[a];
				oldLam[a] = lambda[a];
				oldWeights[a] = weights[a];
			}
			
			iterate();
			
			int counter = 0;
			for (int a = 0; a < mu.length; a++) {
				
				double epsilon = 0.0005;
				if (converges(oldMu[a], mu[a], epsilon)
						&& converges(oldWeights[a], weights[a], epsilon)
						&& converges(oldLam[a], lambda[a], epsilon)) {
					counter++;
				}
			}
			converged = counter == mu.length;
			
			iter++;
			int maxIter = 20;
			if (iter >= maxIter) {
				break;
			}
		}
		
	}
	
	/**
	 * Determine if the EM algorithm has converged
	 *
	 * @param value1  a double representing the value before EM algorithm
	 * @param value2  a double representing the value after EM algorithm
	 * @param epsilon a double representing the maximum difference allowed to be considered converged
	 * @return a boolean indicating whether the values have converged
	 */
	private boolean converges(double value1, double value2, double epsilon) {
		return Math.abs(value1 - value2) <= epsilon;
	}
	
	/**
	 * Access the weighted density of the multi-variate distribution
	 *
	 * @param x      a double representing the read length
	 * @param mean   a double representing the distribution mean
	 * @param lambda  a double representing the distribution standard deviation
	 * @param weight a double representing the distribution weight
	 * @return a double representing the weighted density
	 */
	private double getWeightedDensity(double x, double mean, double lambda, double weight) {
		NormalDistribution dis = new NormalDistribution(mean, lambda);
		return weight * dis.density(x);
	}
	
	/**
	 * Return the index of the largest value in an array of doubles
	 *
	 * @param data an Array of doubles
	 * @return an integer representing thre index of the largest value in the inputted array
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
}
