package Node;

/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
	private String score3 = "";
	private TagNode summit=null;
	private TagNode upstream=null;
	private TagNode downstream=null;
	
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
	 * Set a summit as another TagNode
	 * @param t a TagNode representing the summit
	 */
	public void setSummit(TagNode t){
		summit=t;
	}
	/**
	 * Access the summit 
	 * @return a TagNode representing the summit
	 */
	public TagNode getSummit(){return summit;}
	
	/**
	 * Set the third score variable
	 * @param s a double representing a score
	 */
	public void setScore3(String s){score3 = s;}
	/**
	 * Access third score variable
	 * @return a double representing a score
	 */
	public String getScore3(){return score3;}
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
	 * Access textual description of the entry, including score
	 * @return a String representing textual description of the entry
	 */
	public String toString3(){
		String ans = CHROM+"\t"+BP_START+"\t"+BP_STOP+"\t"+score2;
		return ans;
	}
	/**
	 * Access textual description of the entry, as a scored bedgraph 
	 * @return a String representing the scored bedgraph entry
	 */
	public String toString_ScoredBdg(){
		String ans = CHROM+"\t"+BP_START+"\t"+BP_STOP+"\t"+"E"+(int)score2+"\t"+score3;
		return ans;
	}
	/**
	 * Access textual description of the entry, as a scored summit 
	 * @return a String representing the scored summit entry
	 */
	public String toString_ScoredSummit(){
		
		String ans = CHROM+"\t"+this.getSummit().getStart()+"\t"+this.getSummit().getStop()+"\t"+uniqID+"\t"+score3;
		return ans;
	}
	/**
	 * Access textual description of the entry, as a HMMR gappedPeak 
	 * @return a String representing the scored summit entry
	 */
	public String toString_gappedPeak(){
		
		if (this.upstream == null){
			upstream = this;
		}
		
		if (this.downstream == null){
			downstream=this;
		}
		
		String middleValues = "3"+"\t"+
				"1"+","+getLength()+","+
					"1"+"\t"+"0,"+
				(getStart()-upstream.getStart())+","+
					((downstream.getStop()-upstream.getStart())-1);
		
		String ans = CHROM+"\t"+upstream.getStart()+"\t"+
				downstream.getStop()+"\t"+uniqID+"\t"+".\t."+"\t"+BP_START+"\t"+
				BP_STOP+"\t"+"255,0,0"+"\t"+middleValues+"\t"+score3+"\t"+"-1\t-1";
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
	 * Set an upstream tagnode
	 * @param A TagNode representing the upstream feature
	 */
	public void setUpstream(TagNode u){upstream =u;}
	/**
	 * Get the upstream feature
	 * @return A TagNode representing the upstream feature
	 */
	public TagNode getUpstream(){return upstream;}
	/**
	 * Set an upstream tagnode
	 * @param A TagNode representing the upstream feature
	 */
	public void setDownstream(TagNode d){downstream =d;}
	/**
	 * Get the upstream feature
	 * @return A TagNode representing the upstream feature
	 */
	public TagNode getDownstream(){return downstream;}
	/**
	 * Method for comparing and sorting TagNode
	 */
	public static Comparator<TagNode> basepairComparator = new Comparator<TagNode>() {
        public int compare(TagNode node1, TagNode node2) {
            double PeakScore1 = node1.getStart();
            double PeakScore2 = node2.getStart();
            if (PeakScore1 < PeakScore2) return -1;
            else if(PeakScore1 == PeakScore2) return 0;
            else return 1;
    }
    };
    
    
}
