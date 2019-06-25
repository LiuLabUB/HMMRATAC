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

	private double[] weights;
	private double[] mu;
	private double[] lamda;
	private double[] data;
	private double epsilon = 0.0005;
	private int maxIter=20;
	private double jump = 1.5;

	/**
	 * Constructor for generating new HMMR_EM object
	 * @param w an Array of doubles representing the weights of the distribution
	 * @param m an Array of doubles representing the means of the distribution
	 * @param l an Array of doubles representing the standard deviations of the distribution
	 * @param d an Array of doubles representing the read length data
	 */
	public HMMR_EM(double[] w, double[] m,double[] l,double[] d){
		weights=w;
		mu=m;
		lamda=l;
		data=d;
	}
	/**
	 * Access the means
	 * @return an Array of doubles representing the means
	 */
	public double[] getMeans(){return mu;}
	/**
	 * Access the standard deviations
	 * @return an Array of doubles representing the standard deviations
	 */
	public double[] getLamda(){return lamda;}
	/**
	 * Perform one iteration of the EM algorithm
	 */
	public void iterate(){
		
		double[] temp = new double[mu.length];
		double[] counter = new double[mu.length];
		double total = 0.0;
		Mean[] means = new Mean[mu.length];
		StandardDeviation[] std = new StandardDeviation[mu.length];
		for (int i = 0; i < means.length;i++){
			means[i] = new Mean();
			std[i] = new StandardDeviation();
		}
		for (int i = 0;i < data.length;i++){
			
			for (int l = 0;l < mu.length;l++){
				temp[l] = getWeightedDensity(data[i],mu[l],lamda[l],weights[l]);
				
			}
			int index = returnGreater(temp);
			if (index != -1){
				means[index].increment(data[i]);
				std[index].increment(data[i]);
				counter[index]+=1.0;
				total+=1.0;
			}
		}
		for (int i = 0;i < mu.length;i++){
			mu[i] = mu[i] + (jump *(means[i].getResult()-mu[i]));
			lamda[i] = lamda[i] + (jump *(std[i].getResult() - lamda[i]));
			weights[i] = counter[i] / total;
		}
			
		
		
	}
	/**
	 * Iterate through the EM algorithm until data convergence
	 */
	public void learn(){
		
		double[] oldMu = new double[mu.length];
		double[] oldWeights = new double[mu.length];
		double[] oldLam = new double[mu.length];
		boolean converged = false;
		int iter = 0;
		while(!converged){
			
			for (int a = 0;a < mu.length;a++){
				oldMu[a] = mu[a];
				oldLam[a] = lamda[a];
				oldWeights[a] = weights[a];
			}
			
			
			iterate();
			
			int counter = 0;
			for (int a = 0;a < mu.length;a++){
				
				if(converges(oldMu[a],mu[a],epsilon) && converges(oldWeights[a],weights[a],epsilon)
						&& converges(oldLam[a],lamda[a],epsilon)){
					counter+=1;
				}
			}
			if (counter == mu.length){
				
				converged=true;
			}
			iter+=1;
			if (iter >= maxIter){
				break;
			}
			//Output values during iterations
			/*
			for (int a = 0;a < mu.length;a++){
				System.out.println(mu[a]);
				
				System.out.println(lamda[a]);
				System.out.println(weights[a]);
				System.out.println(iter);
			}
			*/
			//System.out.println(iter);
		}
		
	}
	/**
	 * Determine if the EM algorithm has converged
	 * @param value1 a double representing the value before EM algorithm
	 * @param value2 a double representing the value after EM algorithm
	 * @param epsilon a double representing the maximum difference allowed to be considered converged
	 * @return a boolean indicating whether the values have converged
	 */
	private boolean converges(double value1,double value2, double epsilon){
		if (Math.abs(value1 - value2) <= epsilon){
			return true;
		}
		else{
			return false;
		}
	}
	/**
	 * Access the weighted density of the multi-variate distribution
	 * @param x a double representing the read length
	 * @param mean a double representing the distribution mean
	 * @param lamda	a double representing the distribution standard deviation
	 * @param weight a double representing the distribution weight
	 * @return a double representing the weighted density
	 */
	private double getWeightedDensity(double x, double mean, double lamda,double weight){
		
		NormalDistribution dis = new NormalDistribution(mean,lamda);
		return weight * dis.density(x);
		
	}
	/**
	 * Return the index of the largest value in an array of doubles
	 * @param data an Array of doubles 
	 * @return an integer representing thre index of the largest value in the inputted array
	 */
	private int returnGreater(double[] data){
		int largest = -1;
		double greatest = -1.0;
		for(int i = 0;i < data.length;i++){
			if (data[i] > greatest){
				greatest = data[i];
				largest = i;
			}
		}
		for (int i = 0;i < data.length;i++){
			if (i != largest){
				if(data[i] == greatest){
					largest = -1;
				}
			}
		}
		
		return largest;
	}
		
}
