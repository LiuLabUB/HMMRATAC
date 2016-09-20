package Node;

public class ATACClusterNode extends ATACMatrixNode{

	
	private int _cluster;
	
	public ATACClusterNode(String chrom, int pos, double enrich1,
			double enrich2, double enrich3, int index,int cluster) {
		super(chrom, pos, enrich1, enrich2, enrich3, index);
		_cluster = cluster;
	}
	
	public void setCluster(int cluster){
		_cluster = cluster;
	}
	public int getCluster(){
		return _cluster;
	}

	
	
}
