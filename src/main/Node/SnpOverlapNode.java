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
public class SnpOverlapNode {

	private String snpChr = "";
	private int snpStart = 0;
	private int snpStop = 0;
	private String snpName = "";
	private int snpScore = 0;
	private char snpStrand = '\0';
	private String contigChr = "";
	private int contigStart = 0;
	private int contigStop = 0;
	private String contigName = "";
	private int contigScore = 0;
	private char contigStrand = '\0';
	private int overlap = 0;
	
	public SnpOverlapNode(){
		
	}
	
	public SnpOverlapNode(String snpchr,int snpstart,int snpstop,String snpname,int snpscore,char snpstrand,String contigchr,int contigstart,int contigstop,String contigname,int contigscore,char contigstrand,int over){
		snpChr = snpchr;
		snpStart = snpstart;
		snpStop = snpstop;
		snpName = snpname;
		snpScore = snpscore;
		snpStrand = snpstrand;
		contigChr = contigchr;
		contigStart = contigstart;
		contigStop = contigstop;
		contigName = contigname;
		contigScore = contigscore;
		contigStrand = contigstrand;
		overlap = over;
	}
	
	public String getSnpChr(){
		return snpChr;
	}
	public void setSnpChr(String snpchr){
		snpChr = snpchr;
	}
	
	public int getSnpStart(){
		return snpStart;
	}
	public void setSnpStart(int snpstart){
		snpStart = snpstart;
	}
	
	public int getSnpStop(){
		return snpStart;
	}
	public void setSnpStop(int snpstop){
		snpStop = snpstop;
	}
	
	public String getSnpName(){
		return snpName;
	}
	public void setSnpName(String snpname){
		snpName= snpname;
	}
	
	public char getSnpStrand(){
		return snpStrand;
	}
	public void setSnpStrand(char snpstrand){
		snpStrand = snpstrand;
	}
	
	public String getContigChr(){
		return contigChr;
	}
	public void setContigChr(String contigchr){
		contigChr = contigchr;
	}
	
	public int getContigStart(){
		return contigStart;
	}
	public void setContigStart(int contigstart){
		contigStart=contigstart;
	}
	
	public int getContigStop(){
		return contigStop;
	}
	public void setContigStop(int contigstop){
		contigStop = contigstop;
	}
	
	public String getContigName(){
		return contigName;
	}
	public void setContigName(String contigname){
		contigName = contigname;
	}
	
	public char getContigStrand(){
		return contigStrand;
	}
	public void setContigStrand(char contigstrand){
		contigStrand = contigstrand;
	}
	
	public int getOverlap(){
		return overlap;
	}
	public void setOverlap(int over){
		overlap=over;
	}
}
