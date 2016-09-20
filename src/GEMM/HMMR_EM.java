package GEMM;


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

	
	public HMMR_EM(double[] w, double[] m,double[] l,double[] d){
		weights=w;
		mu=m;
		lamda=l;
		data=d;
	}
	public double[] getMeans(){return mu;}
	public double[] getLamda(){return lamda;}
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
			System.out.println(iter);
		}
		
	}
	private boolean converges(double value1,double value2, double epsilon){
		if (Math.abs(value1 - value2) <= epsilon){
			return true;
		}
		else{
			return false;
		}
	}
	private double getWeightedDensity(double x, double mean, double lamda,double weight){
		
		NormalDistribution dis = new NormalDistribution(mean,lamda);
		return weight * dis.density(x);
		
	}
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
