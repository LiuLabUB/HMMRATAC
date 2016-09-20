package Node;
import java.util.Vector;


public class MyTagNode {

	private String Chrom = "";
	private int Start = 0;
	private int Stop = 0;
	private Vector<TagNode> Data = new Vector<TagNode>();
	
	
	public MyTagNode(){
		
	}
	public MyTagNode(String chr,int start,int stop,Vector<TagNode> data){
		Start = start;
		Stop = stop;
		Chrom = chr;
		Data = data;
		
	}
	public void setTemp(Vector<TagNode> data){
		Data = data;
	}
	
	public Vector<TagNode> getData(){
		return Data;
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
	
	public String getChrom(){
		return Chrom;
	}
	public void setChrom(String chr){
		Chrom = chr;
	}
}
