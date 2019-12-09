
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

public class BedGraphNode extends Node.TagNode{
	
	private String _chrom;
	private int _start;
	private int _stop;
	private int _score;
	private double _score2;
	private String _score3;
	
	public BedGraphNode(String chr,int start,int stop,int score){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score = score;
	}
	public BedGraphNode(String chr,int start,int stop,double score){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score2 = score;
	}
	public BedGraphNode(){
		
	}
    public BedGraphNode(String chr,int start,int stop,double score,String s){
		_chrom = chr;
		_start = start;
		_stop = stop;
		_score2 = score;
		_score3 = s;
	}
	
	public String getChrom(){return _chrom;}
	public void setChrom(String chr){_chrom = chr;}
	public int getStart(){return _start;}
	public void setStart(int s){_start = s;}
	public int getStop(){return _stop;}
	public void setStop(int s){_stop = s;}
	public int getScore(){return _score;}
	public void setScore(int s){_score = s;}
	public double getScore2(){return _score2;}
	public void setScore2(double s){_score2 = s;}
	public String getScore3(){return _score3;}
	public void setScore3(String score){_score3 = score;}
}
