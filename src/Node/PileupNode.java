package Node;

public class PileupNode {
	
	//private String _chrom;
	private int _base;
	private double _score;
	private int _score2;
	public PileupNode(int base, double score){
	//	_chrom = chr;
		_base = base;
		_score = score;
	}
	public PileupNode(){
		
	}
	public PileupNode(int base,int score){
		_base = base;
		_score2 = score;
	}
	/*
	public String getChrom(){
		return _chrom;
	}
	public void setChrom(String chr){
		_chrom = chr;
	}
	*/
	public int getBase(){
		return _base;
	}
	public void setBase(int base){
		_base = base;
	}
	
	public double getScore(){
		return _score;
	}
	public void setScore(double score){
		_score = score;
	}
	public int getScore2(){
		return _score2;
	}
	public void setScore2(int score){
		_score2 = score;
	}
}
