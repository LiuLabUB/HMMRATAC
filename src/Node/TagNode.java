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
	
	public TagNode () {
		
	}

	public TagNode (String chr, int start, int stop,String uniqid,int score,char str) {
		CHROM = chr;
		BP_START = start;
		BP_STOP = stop;
		uniqID = uniqid;
		strand = str;
		Score = score;
		
	}
	public TagNode(String chr, int start, int stop) {
		
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
	}
	public TagNode(String chr,int start,int stop,int score){
		CHROM = chr;
		BP_START = start;
		BP_STOP = stop;
		Score = score;
		
	}
	public TagNode(String chr,int start,int stop,double score){
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
		score2 = score;
	}
	public void setScore3(double s){score3 = s;}
	public double getScore3(){return score3;}
	public void setScore2(double score){
		score2=score;
	}
	public double getScore2(){
		return score2;
	}
	public int getLength(){
		return BP_STOP - BP_START;
	}

	public char getcharStrand(){
		return strand;
	}
	public void setcharStrand(char str){
		strand = str;
	}
	
	public int getOverlap() {
		return overlap;
	}
	
	public void setOverlap(int newoverlap) {
		overlap = newoverlap;
	}
	
	public SAMRecord getRecord() {
		return record;
	}
	public int getScore(){
		return Score;
	}
	public void setScore(int score){
		Score = score;
	}
	
	public void setRecord(SAMRecord newrecord) {
		record = newrecord;
	}
	
	public String getID() {
		return uniqID;		
	}
	
	public void setID(String newid) {
		uniqID = newid;		
	}
	
	public String getChrom() {
		return CHROM;		
	}
	
	public void setChrom(String newchr) {
		CHROM = newchr;		
	}
		
	public int getStop() {
		return BP_STOP;		
	}
	
	public void setStop(int newbp) {
		BP_STOP = newbp;		
	}
	
	public int getStart() {
		return BP_START;		
	}
	
	public void setStart(int newbp) {
		BP_START = newbp;		
	}
	
	public boolean getStrand() {
		return Fstrand;		
	}
	
	public void setStrand(boolean newstrand) {
		Fstrand = newstrand;		
	}
	
	public static Comparator<TagNode> basepairComparator = new Comparator<TagNode>() {
        public int compare(TagNode node1, TagNode node2) {
            double PeakScore1 = node1.getStart();
            double PeakScore2 = node2.getStart();
            if (PeakScore1 < PeakScore2) return -1;
            else return 1;
    }
    };
    
    
}