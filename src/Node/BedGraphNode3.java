package Node;

public class BedGraphNode3 extends BEDGraphWithIndexNode{

	
	private int _cluster;
	public BedGraphNode3(String chr, int start, int stop, double score,
			int index,int cluster) {
		super(chr, start, stop, score, index);
		_cluster = cluster;
		
	}

	public BedGraphNode3(String chr, int start, int stop, int score,
			int index,int cluster) {
		super(chr, start, stop, score, index);
		_cluster = cluster;
	}
	public BedGraphNode3(){
		
	}
	public int getCluster(){return _cluster;}
	public void setCluster(int c){_cluster = c;}
	
	
}
