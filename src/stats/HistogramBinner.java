package stats;
public class HistogramBinner {
	
	private double BinStart = 0.0;
	private double BinStop = 0.0;
	private int BinNumber = 0;
	
	public HistogramBinner(){
		
	}

	public HistogramBinner(double binstart,double binstop,int binnumber){
		BinStart = binstart;
		BinStop = binstop;
		BinNumber = binnumber;
	}
	
	public double getStart(){
		return BinStart;
	}
	public void setStart(double binstart){
		BinStart = binstart;
	}
	
	public double getStop(){
		return BinStop;
	}
	public void setStop(double binstop){
		BinStop = binstop;
	}
	
	public int getNumber(){
		return BinNumber;
	}
	public void setNumber(int binnumber){
		BinNumber = binnumber;
	}
}
