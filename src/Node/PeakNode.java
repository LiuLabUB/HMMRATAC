package Node;
public class PeakNode {

	private int Summit = 0;
	private int Start = 0;
	private int Stop = 0;
	
	//public PeakNode(){
		
	//}
	
	public PeakNode(int summit,int start,int stop){
		Start = start;
		Summit = summit;
		Stop = stop;
	}
	
	public int getSummit(){
		return Summit;
	}
	public void setSummit(int summit){
		Summit = summit;
	}
	
	public int getStart(){
		return Start;
	}
	public void setStart(int start){
		Start = start;
	}
	
	public int getStop(){
		return Stop;
	}
	public void setStop(int stop){
		Stop = stop;
	}
}
