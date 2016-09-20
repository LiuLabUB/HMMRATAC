package Node;

public class PileupNode2 extends PileupNode{

	
	private String chrom;
	private int index;
	private int moreThanMinReads;
	private int partOfShortPeak;
	private int partOfMonoPeak;
	private double probability;
	
	public PileupNode2(int base, double score,String chr,int i) {
		super(base, score);
		chrom = chr;
		index = i;
	}
	public PileupNode2(){
		
	}
	public PileupNode2(int base,int score,String chr){
		super(base,score);
		chrom = chr;
	}
	public PileupNode2(int base,int score,String chr,double prob){
		super(base,score);
		chrom = chr;
		probability = prob;
	}
	public PileupNode2(int base,double score,String chr){
		super(base,score);
		chrom = chr;
	}
	public PileupNode2(int base, double score,String chr,int i,int more,int partSh, int partMo) {
		super(base, score);
		chrom = chr;
		index = i;
		moreThanMinReads = more;
		partOfShortPeak = partSh;
		partOfMonoPeak = partMo;
	}
	
	public void setIndex(int i){
		index = i;
	}
	public int getIndex(){
		return index;
	}
	public void setChrom(String chr){
		chrom = chr;
	}
	public String getChrom(){
		return chrom;
	}
	public void setMoreThanMin(int more){
		moreThanMinReads = more; 
	}
	public int getMoreThanMin(){
		return moreThanMinReads;
	}
	public void setPartOfShort(int part){
		partOfShortPeak = part;
	}
	public int getPartOfShort(){
		return partOfShortPeak;
	}
	public void setPartOfMono(int part){
		partOfMonoPeak= part;
	}
	public int getPartOfMono(){
		return partOfMonoPeak;
	}
	public double getProb(){return probability;}
	public void setProb(double prob){probability = prob;}
}
