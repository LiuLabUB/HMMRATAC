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

public class AlleleNode {
	private String SnpName = "";
	private String ConName = "";
	private int Position = 0;
	private char NT = '\0';
	private int ConStart = 0;
	private int ConStop = 0;
	private String Header = "";
	private StringBuilder Seq;
	private String Quality = "";
	
	public AlleleNode(){
		
	}
	
	public AlleleNode(String snp,String con,int pos,char nt,int start,int stop,String head,StringBuilder seq, String qual){
		SnpName = snp;
		ConName = con;
		Position = pos;
		NT = nt;
		ConStart = start;
		ConStop = stop;
		Header = head;
		Seq = seq;
		Quality = qual;
	}
	public int getPos(){
		return Position;
	}
	public void setPos(int pos){
		Position = pos;
	}
	public int getStart(){
		return ConStart;
	}
	public void setStart(int start){
		ConStart = start;
	}
	public int getStop(){
		return ConStop;
	}
	public void setStop(int stop){
		ConStop = stop;
	}
	
	public String getSnpName(){
		return SnpName;
	}
	public void setSnp(String snp){
		SnpName = snp;
	}
	
	public String getCon(){
		return ConName;
	}
	public void setCon(String con){
		ConName = con;
	}
	
	public char getNT(){
		return NT;
	}
	public void setNT(char nt){
		NT = nt;
	}
	
	public String getHeader(){
		return Header;
	}
	public void setHeader(String head){
		Header = head;
	}
	
	public StringBuilder getSeq(){
		return Seq;
	}
	public void setSeq(StringBuilder seq){
		Seq = seq;
	}
	
	public String getQuality(){
		return Quality;
	}
	public void setQual(String qual){
		Quality = qual;
	}

}
