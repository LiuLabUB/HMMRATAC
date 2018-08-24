package Node;

public class PileupNode {
	
	
	private int _base;
	private double _score;
	private int _score2;
	/**
	 * Constructor for creating new PileupNode
	 * @param base an integer representing the base position
	 * @param score a double representing the base's score
	 */
	public PileupNode(int base, double score){
	//	_chrom = chr;
		_base = base;
		_score = score;
	}
	/**
	 * Default Constructor
	 */
	public PileupNode(){
		
	}
	/**
	 * Constructor for creating new PileupNode
	 * @param base an integer representing the base position
	 * @param score a integer representing the base's score
	 */
	public PileupNode(int base,int score){
		_base = base;
		_score2 = score;
	}
	/**
	 * Access the base position
	 * @return an integer representing the base position
	 */
	public int getBase(){
		return _base;
	}
	/**
	 * Set the base position
	 * @param base an integer representing the base position
	 */
	public void setBase(int base){
		_base = base;
	}
	/**
	 * Access the score
	 * @return a double representing the score at the position
	 */
	public double getScore(){
		return _score;
	}
	/**
	 * Set the score
	 * @param score a double representing the positions score
	 */
	public void setScore(double score){
		_score = score;
	}
	/**
	 * Access the score
	 * @return an integer representing the score at the position
	 */
	public int getScore2(){
		return _score2;
	}
	/**
	 * Set the score
	 * @param score an integer representing the positions score
	 */
	public void setScore2(int score){
		_score2 = score;
	}
}
