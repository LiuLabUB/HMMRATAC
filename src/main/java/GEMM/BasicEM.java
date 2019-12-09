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

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import stats.RectangularSmooth;

public class BasicEM {

	private double[] mode;
	private double[] weights;
	private double[] mu;
	private double[] lamda;
	private double[] data;
	private int maxIter = 10;
	private double epsilon = 0.0005;
	private int numIter = 0;
	private double jump = 1.5;
	
	public BasicEM(double[] m,double[] w,double[] means,double[] l,double[] d){
		mode = m;
		weights = w;
		mu = means;
		lamda = l;
		data = d;
	}
	public void setMaxIter(int m){maxIter = m;}
	public void setEpsilon(double e){epsilon = e;}
	
	public int getNumberIterations(){return numIter;}
	public double[] getMeans(){return mu;}
	public double[] getWeights(){return weights;}
	public double[] getLamda(){return lamda;}
	
	public void learn(){
		
		int i;
		double[] oldMu = new double[mode.length];
		double[] oldWeights = new double[mode.length];
		double[] oldLam = new double[mode.length];
		
		for ( i = 0; i < maxIter;i++){
			
			for (int a = 0;a < mode.length;a++){
				oldMu[a] = mu[a];
				oldLam[a] = lamda[a];
				oldWeights[a] = weights[a];
			}
			
			
			iterate();
			
			//Out statement 
			//for(int x = 0;x < mode.length;x++){
				//System.out.println(weights[x]+"\t"+mu[x]+"\t"+lamda[x]);
			//}
			
			int counter = 0;
			for (int a = 0;a < mode.length;a++){
				//System.out.println(oldMu[a]+"\t"+mu[a]);
				if(converges(oldMu[a],mu[a],epsilon) && converges(oldWeights[a],weights[a],epsilon)
						&& converges(oldLam[a],lamda[a],epsilon)){
					counter+=1;
				}
			}
			if (counter == mode.length){
				
				break;
			}
			
		}
		numIter = i;
	}
	
	private boolean converges(double value1,double value2, double epsilon){
		if (Math.abs(value1 - value2) <= epsilon){
			return true;
		}
		else{
			return false;
		}
	}
	
	public void iterate(){
		ArrayList<ArrayList<Double>> separatedData = new ArrayList<ArrayList<Double>>();
		for (int i = 0;i < mode.length;i++)
			separatedData.add(new ArrayList<Double>());
		double[] temp = new double[mode.length];
		double counter = 0.0;
		for (int i = 0;i < data.length;i++){
			
			for (int l = 0;l < mode.length;l++){
				temp[l] = getWeightedDensity(mode[l],data[i],mu[l],lamda[l],weights[l]);
				
			}
			int index = returnGreater(temp);
			if (index != -1){
				separatedData.get(index).add(data[i]);
				counter += 1.0;
			}
		}
		
		for (int l = 0;l < mode.length;l++){
			double[] d = convertFromList(separatedData.get(l));
			if (mode[l] == 1){
				double oldMu1 = mu[l];
				double oldLamda1 = lamda[l];
				mu[l] = getArithmicMedian(d);
				lamda[l] = updateLaplaceBeta(d);
				
			}
			else if(mode[l] == 2){
				//double oldMu2 = mu[l];
				//double oldLamda2 = lamda[l];
				mu[l] = mu[l] + (jump*( getArithmicMean(d)-mu[l]) );
				lamda[l] = lamda[l] + (jump*(getArithmicStdDev(d)-lamda[l]));
				//double muDiff = Math.abs(mu[l] - oldMu2);
				//double lamdaDiff = Math.abs(lamda[l] - oldLamda2);
				//mu[l] = mu[l] * (muDiff * jump);
				//lamda[l] = lamda[l] * (lamdaDiff * jump);
			}
			else if(mode[l] == 0.5){
				double oldMu3 = mu[l];
				double oldLamda3 = lamda[l];
				mu[l] = getArithmicMean(d);
				lamda[l] = lamda[l];
			}
			else if(mode[l] == 0){
				double oldMu4 = mu[l];
				double oldLamda4 = lamda[l];
				mu[l] = mu[l];
				lamda[l] = updateLaplaceBeta(d);
			}
			else if(mode[l] == 3){
				double oldMu5 = mu[l];
				double oldLamda5 = lamda[l];
				mu[l] = mu[l];
				lamda[l] = lamda[l];
			}
			double xml = (double) d.length / counter;
			weights[l] = weights[l] + (jump * (xml-weights[l]));
		}
		
	}
	/*
	 * returns the index of data that is the greatest value or -1 if the largest value in data is duplicated
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
	private double getWeightedDensity(double p, double x, double mean, double lamda,double weight){
		double density = 0.0;
		
		 if (p == 2){
			NormalDistribution dis = new NormalDistribution(mean,lamda);
			density = weight * dis.density(x);
		}
		
		else if (p == 0.5 || p == 3){
			ExponentialDistribution dis = new ExponentialDistribution(mean);
			density = weight * dis.density(x);
		}
		
		return density;
	}
	private double getDensity(double p,double x,double mean,double lamda){
		double density = 0.0;
		
		if (p == 2){
			NormalDistribution dis = new NormalDistribution(mean,lamda);
			density = dis.density(x);
		}
		
		if (p == 0.5){
			ExponentialDistribution dis = new ExponentialDistribution(mean);
			density = dis.density(x);
		}
		
		return density;
	}
	private double[] convertFromList(ArrayList<Double> data){
		double[] newData = new double[data.size()];
		for (int i = 0;i < data.size();i++){
			newData[i] = data.get(i);
		}
		return newData;
	}
	private double getArithmicMean(double[] data){
		return new Mean().evaluate(data);
	}
	private double getArithmicStdDev(double[] data){
		return new StandardDeviation().evaluate(data);
	}
	private double getArithmicMedian(double[] data){
		return new Percentile().evaluate(data, 50.0);
	}
	private double updateLaplaceBeta(double[] data){
		double mean = getArithmicMedian(data);
		double sum = 0.0;
		double counter = 0.0;
		for (int i = 0; i < data.length;i++){
			sum += Math.abs(data[i] - mean);
			counter += 1.0;
		}
		return sum/counter;
	}
	
	public static void main(String[] args){
		String d = args[0];
		Scanner inFile = null;
		try{
			inFile = new Scanner((Readable) new FileReader(d));
		}
		catch(FileNotFoundException e){
			e.printStackTrace();
		}
		int counter = 0;
		ArrayList<Double> temp = new ArrayList<Double>();
		while (inFile.hasNextLine()){
			String line = inFile.nextLine();
			String[] feat = line.split("//s+");
			double value = Double.parseDouble(feat[0]);
			temp.add(value);
		}
		double[] data = new double[temp.size()];
		double[] hist = new double[1001];
		for (int i = 0; i < temp.size();i++){
			data[i] = temp.get(i);
			//System.out.println(temp.get(i).toString().split("\\.").length);
			hist[Integer.parseInt(temp.get(i).toString().split("\\.")[0])] += 1;
			
		}
		temp = null;
		/*for (double value : hist){
			System.out.println(value);
		}*/
		
		//idea is to use the local max and local mins to identify the initial parameters. short initial mean will be readlength
		//mono initial could be first local max (certain distance from short) and di/tri could be 2x or 3x of mono etc
		/*
		RectangularSmooth smoother = new RectangularSmooth(hist,30);
		hist = smoother.getsmoothedData();
		
		int[] max = getLocalExtrema(hist,"greater",1);
		int[] min = getLocalExtrema(hist,"less than",1);
		for (int m : max)
			System.out.println("Max\t"+m);
		for(int m : min)
			System.out.println("Min\t"+m);
		*/
		
		double[] weights = new double[3];
		double[] means = new double[3];
		double[] lam = new double[3];
		double[] mode = new double[3];
		
		/*
		means[0] = 50;
		lam[0] = 5;
		mode[0] = 0;
		weights[0] = 0.25;
		*/
		
		means[0] = 150; means[1] = 300; means[2] = 450;
		lam[0] = lam[1] = lam[2] = 20;
		mode[0] = mode[1] = mode[2] = 2;
		//weights[1] = 0.41; weights[2] = 0.17; weights[3] = 0.17;
		
		
		double cutone = (means[1] - means[0])/2 + means[0];
		double cuttwo = (means[2] - means[1])/2 + means[1];
		int counter1=0;
		int counter2=0;
		int counter3=0;
		for (int i = 0;i < data.length;i++){
			if (data[i] < cutone)
				counter1++;
			else if(data[i] >= cutone && data[i] <= cuttwo)
				counter2++;
			else
				counter3++;
		}
		weights[0] = (double)counter1/(double)data.length; 
		weights[1] = (double)counter2/(double)data.length;
		weights[2] = (double)counter3/(double)data.length;
		
		for(int i = 0;i < mode.length;i++){
			System.out.println(weights[i]+"\t"+means[i]+"\t"+lam[i]);
		}
		
		BasicEM em = new BasicEM(mode,weights,means,lam,data);
		em.setMaxIter(50);
		em.learn();
		/*
		for (int a = 0;a < 10;a++){
			em.iterate();
			//System.out.println(em.getNumberIterations());
			means = em.getMeans();
			weights = em.getWeights();
			lam = em.getLamda();
			for(int i = 0;i < mode.length;i++){
				System.out.println(weights[i]+"\t"+means[i]+"\t"+lam[i]);
			}
		}*/
		
	}
	public static int[] getInflectionPoints(double[] smoothed,int order){
		int[] local;
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i = order; i < smoothed.length - order;i++){
			double upSlope = Math.abs(smoothed[i-order] - smoothed[i]) / (double)order+1.0;
			double downSlope = Math.abs(smoothed[i] - smoothed[i+order]) / (double)order+1.0;
			if (upSlope != downSlope){
				temp.add(i);
			}
		}
		local = new int[temp.size()];
		for (int i = 0;i < temp.size();i++)
			local[i] = temp.get(i);
		return local;
	}
	public static int[] getLocalExtrema(double[] smoothed,String comparator,int order){
		int[] local;// = null;
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for (int i = (order);i < smoothed.length - (order);i++){
			if (comparator.toLowerCase().contains("great")){
				if (Math.max(smoothed[i], smoothed[i-order]) == smoothed[i] &&
						Math.max(smoothed[i], smoothed[i+order]) == smoothed[i]&&
						smoothed[i] != smoothed[i+order]&&
						smoothed[i] != smoothed[i-order]){
					temp.add(i);
				}
			}
			else if(comparator.toLowerCase().contains("less")){
				if (Math.min(smoothed[i], smoothed[i-order]) == smoothed[i] &&
						Math.min(smoothed[i], smoothed[i+order]) == smoothed[i]&&
								smoothed[i] != smoothed[i+order]&&
								smoothed[i] != smoothed[i-order]){
					temp.add(i);
				}
			}
		}
		local = new int[temp.size()];
		for (int i = 0; i < temp.size();i++){
			local[i] = temp.get(i);
		}
		
		return local;
	}
	public static int[] getLocalMaxMin(double[] _smoothed){
		int[] local = new int[_smoothed.length];
		for (int i =2;i < _smoothed.length-2;i++){
			//Calculates all local maxima and local minima points in the smoothed pileups
			
			double pointone = _smoothed[i-2];
			double pointtwo = _smoothed[i-1];
			double pointthree = _smoothed[i];
			double pointfour = _smoothed[i+1];
			double pointfive = _smoothed[i+2];
			
			double slopepointtwo = (pointthree-pointone)/((i)-(i-2));
			double slopepointthree = (pointfour - pointtwo)/((i+1)-(i-1));
			double slopepointfour = (pointfive-pointthree)/((i+2)-(i));
			
			if (slopepointthree == 0){
				double secondslopepointthree = (slopepointfour - slopepointtwo)/((i+1)-(i-1));
				if (secondslopepointthree > 0){local[i]--;}
				if (secondslopepointthree < 0){local[i]++;}
			}
		}
		return local;
	}
}
