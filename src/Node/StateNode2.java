package Node;

public class StateNode2 extends StateNode{

	
	private double emissionProb;
	private double transitionProb;
	private double zscore;
	public StateNode2(String c, int p, int s,double ep,double tp) {
		super(c, p, s);
		emissionProb = ep;
		transitionProb = tp;
	}
	public StateNode2(String c,int p,int s,double ep,double tp,double z){
		super(c, p, s);
		emissionProb = ep;
		transitionProb = tp;
		zscore = z;
	}
	
	public void setEP(double ep){
		emissionProb = ep;
	}
	public double getEP(){
		return emissionProb;
	}
	public void setTP(double tp){
		transitionProb = tp;
	}
	public double getTP(){
		return transitionProb;
	}
	public void setZScore(double z){zscore = z;}
	public double getZScore(){return zscore;}
}
