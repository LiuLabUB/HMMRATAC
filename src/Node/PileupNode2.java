package Node;

public class PileupNode2 extends PileupNode{

	
	private String chrom;
	
	/**
	 * Constructor for creating PileupNode2 object
	 * @param base an integer representing the base position
	 * @param score an integer representing the position's score
	 * @param chr a String representing the chromosome name
	 */
	public PileupNode2(int base,int score,String chr){
		super(base,score);
		chrom = chr;
	}
	/**
	 * Constructor for creating PileupNode2 object
	 * @param base an integer representing the base position
	 * @param score a double representing the position's score
	 * @param chr a String representing the chromosome name
	 */
	public PileupNode2(int base,double score,String chr){
		super(base,score);
		chrom = chr;
	}
	
	
	/**
	 * Set the chromsome name
	 * @param chr a String representing the chromosome name
	 */
	public void setChrom(String chr){
		chrom = chr;
	}
	/**
	 * Access the chromosome name
	 * @return a String representing the chromosome na,e
	 */
	public String getChrom(){
		return chrom;
	}
	
}
