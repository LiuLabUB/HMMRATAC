package Node;

import java.util.ArrayList;

public class ScoreNode {
	private double ave;
	private double median;
	private double max;
	private ArrayList<Double> val;
	
	public ScoreNode(double a,double med, double m,ArrayList<Double> v){
		ave=a;median=med;max=m;val=v;
	}
	public double getMean(){return ave;}
	public double getMax(){return max;}
	public double getMedian(){return median;}
	public ArrayList<Double> getvalues(){return val;}
}
