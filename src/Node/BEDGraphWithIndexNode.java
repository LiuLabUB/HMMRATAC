package Node;

public class BEDGraphWithIndexNode extends BedGraphNode{

	
	private int _index;
	public BEDGraphWithIndexNode(String chr, int start, int stop, double score,int index) {
		super(chr, start, stop, score);
		_index = index;
	}

	public BEDGraphWithIndexNode(String chr, int start, int stop, int score,int index) {
		super(chr, start, stop, score);
		_index = index;
	}
	public BEDGraphWithIndexNode(){
		
	}
	
	public int getIndex(){return _index;}
	public void setIndex(int i){_index = i;}
	
}
