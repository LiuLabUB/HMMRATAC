
package Node;

public class BedGraphNode extends Node.TagNode{
	
	private String _chrom;
	private int _start;
	private int _stop;
	private int _score;
	private double _score2;
	private double _score3;
	
	public BedGraphNode(String chr,int start,int stop,int score){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score = score;
	}
	public BedGraphNode(String chr,int start,int stop,double score){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score2 = score;
	}
	public BedGraphNode(){
		
	}
	public BedGraphNode(String chr,int start,int stop,double score,double s){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score2 = score;
		_score3 = s;
	}
	
	public String getChrom(){return _chrom;}
	public void setChrom(String chr){_chrom = chr;}
	public int getStart(){return _start;}
	public void setStart(int s){_start = s;}
	public int getStop(){return _stop;}
	public void setStop(int s){_stop = s;}
	public int getScore(){return _score;}
	public void setScore(int s){_score = s;}
	public double getScore2(){return _score2;}
	public void setScore2(double s){_score2 = s;}
	public double getScore3(){return _score3;}
	public void setScore3(double score){_score3 = score;}
}
