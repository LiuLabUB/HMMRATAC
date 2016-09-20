package Node;

public class DENode extends TagNode{

	private double average1;
	private double average2;
	private double pvalue;
	
	public DENode(String chr,int start, int stop,double ave1,double ave2,double p){
		super(chr,start,stop);
		average1 = ave1;
		average2 = ave2;
		pvalue = p;
	}
	public void setAve1(double ave1){
		average1 = ave1;
	}
	public double getAve1(){
		return average1;
	}
	public void setAve2(double ave2){
		average2 = ave2;
		
	}
	public double getAve2(){
		return average2;
	}
	public void setPValue(double p){
		pvalue = p;
	}
	public double getPValue(){
		return pvalue;
	}
	
}
