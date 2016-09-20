package Node;

public class MatrixNodeForKMeans extends ATACMatrixNode{

	private int dataID;
	
	public MatrixNodeForKMeans(String chrom, int pos, double enrich1,
			double enrich2, double enrich3, int index, int ID) {
		super(chrom, pos, enrich1, enrich2, enrich3, index);
		dataID = ID;
	}
	
	public void setID(int ID){
		dataID = ID;
	}
	public int getID(){
		return dataID;
	}

}
