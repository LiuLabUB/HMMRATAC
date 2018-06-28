package Node;
import java.util.Comparator;

import net.sf.samtools.SAMRecord;

public class TagNode {
	private String uniqID = "";
	private String CHROM = "";
	private int BP_START = 0;
	private int BP_STOP = 0;
	private boolean Fstrand = true;
	private char strand = '\0';
	private int Score = 0;
	private double score2 = 0.0;
	private double score3 = 0.0;
	
	private int overlap = 0;
	
	private SAMRecord record = null;
	
	/**
	 * Default constructor
	 */
	public TagNode () {
		
	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 * @param uniqid A String representing a unique ID for the entry
	 * @param score an integer representing the score for the entry
	 * @param str a Char representing the strand of the entry
	 */
	public TagNode (String chr, int start, int stop,String uniqid,int score,char str) {
		CHROM = chr;
		BP_START = start;
		BP_STOP = stop;
		uniqID = uniqid;
		strand = str;
		Score = score;
		
	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 */
	public TagNode(String chr, int start, int stop) {
		
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 * @param score an integer representing the score for the entry
	 */
	 
	public TagNode(String chr,int start,int stop,int score){
		CHROM = chr;
		BP_START = start;
		BP_STOP = stop;
		Score = score;
		
	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 * @param score a double representing the score for the entry
	 */
	public TagNode(String chr,int start,int stop,double score){
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
		score2 = score;
	}
	/**
	 * Set the third score variable
	 * @param s a double representing a score
	 */
	public void setScore3(double s){score3 = s;}
	/**
	 * Access third score variable
	 * @return a double representing a score
	 */
	public double getScore3(){return score3;}
	/**
	 * Set Score
	 * @param score a double representing the entry's score
	 */
	public void setScore2(double score){
		score2=score;
	}
	/**
	 * Access the score
	 * @return a double representing the entry's score
	 */
	public double getScore2(){
		return score2;
	}
	/**
	 * Access the length of the entry
	 * @return an integer representing the length of the entry
	 */
	public int getLength(){
		return BP_STOP - BP_START;
	}
	/**
	 * Access textual description of the entry
	 * @return a String representing textual description of the entry
	 */
	public String toString(){
		String ans = CHROM+"\t"+BP_START+"\t"+BP_STOP;
		return ans;
	}
	/**
	 * Access textual description of the entry
	 * @return a String representing textual description of the entry
	 */
	public String toString2(){
		String ans = CHROM+"\t"+BP_START+"\t"+BP_STOP+"\t"+"E"+(int)score2;
		return ans;
	}
	/**
	 * Access the entry's strand
	 * @return a char representing the strand for the entry
	 */
	public char getcharStrand(){
		return strand;
	}
	/**
	 * Set the entry's strand
	 * @param str a char representing the entry's strand
	 */
	public void setcharStrand(char str){
		strand = str;
	}
	/**
	 * Access the overlap
	 * @return an integer representing the overlap
	 */
	public int getOverlap() {
		return overlap;
	}
	/**
	 * Set the overlap
	 * @param newoverlap an integer representing the overlap
	 */
	public void setOverlap(int newoverlap) {
		overlap = newoverlap;
	}
	/**
	 * Access the SAM record
	 * @return a SAMRecord for the entry
	 */
	public SAMRecord getRecord() {
		return record;
	}
	/**
	 * Access the integer score
	 * @return an integer representing the entry's score
	 */
	public int getScore(){
		return Score;
	}
	/**
	 * Set the integer score
	 * @param score an integer representing the entry's score
	 */
	public void setScore(int score){
		Score = score;
	}
	/**
	 * Set the SAM record
	 * @param newrecord a SAMRecord for the entry
	 */
	public void setRecord(SAMRecord newrecord) {
		record = newrecord;
	}
	/**
	 * Access the unique ID
	 * @return a String representing the entry's unique ID
	 */
	public String getID() {
		return uniqID;		
	}
	/**
	 * Set the unique ID
	 * @param newid a String representing the entry's unique ID
	 */
	public void setID(String newid) {
		uniqID = newid;		
	}
	/**
	 * Access the chromosome name
	 * @return a String representing the chromosome name
	 */
	public String getChrom() {
		return CHROM;		
	}
	/**
	 * Set the chromosome name
	 * @param newchr a String representing the chromosome name
	 */
	public void setChrom(String newchr) {
		CHROM = newchr;		
	}
	/**
	 * Access the stop position	
	 * @return an integer representing the stop position
	 */
	public int getStop() {
		return BP_STOP;		
	}
	/**
	 * Set the stop position
	 * @param newbp an integer representing the stop position
	 */
	public void setStop(int newbp) {
		BP_STOP = newbp;		
	}
	/**
	 * Access the start position
	 * @return an integer representing the start position
	 */
	public int getStart() {
		return BP_START;		
	}
	/**
	 * Set the start position
	 * @param newbp an integer representing the start position
	 */
	public void setStart(int newbp) {
		BP_START = newbp;		
	}
	/**
	 * Access whether the entry is on the forward strand
	 * @return a a boolean to determine if the entry is on the forward strand
	 */
	public boolean getStrand() {
		return Fstrand;		
	}
	/**
	 * Set whether the entry is on the forward strand
	 * @param newstrand a boolean to determine id the entry is on the forward strand
	 */
	public void setStrand(boolean newstrand) {
		Fstrand = newstrand;		
	}
	/**
	 * Method for comparing and sorting TagNode
	 */
	public static Comparator<TagNode> basepairComparator = new Comparator<TagNode>() {
        public int compare(TagNode node1, TagNode node2) {
            double PeakScore1 = node1.getStart();
            double PeakScore2 = node2.getStart();
            if (PeakScore1 < PeakScore2) return -1;
            else return 1;
    }
    };
    
    
}