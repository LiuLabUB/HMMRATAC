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
public class ThreeSiteBedNode {

	
	//Can handle two or three site bed lines
	private String chrom1;
	private int start1;
	private int stop1;
	private String chrom2;
	private int start2;
	private int stop2;
	private String chrom3;
	private int start3;
	private int stop3;
	private boolean containsThree;
	
	public ThreeSiteBedNode(String c1, int s1, int st1, String c2, int s2, int st2, String c3, int s3, int st3){
		chrom1 = c1;
		start1 = s1;
		stop1 = st1;
		chrom2 = c2;
		start2 = s2;
		stop2 = st2;
		chrom3 = c3;
		start3 = s3;
		stop3 = st3;
		containsThree = true;
	}
	public ThreeSiteBedNode(String c1, int s1, int st1, String c2, int s2, int st2){
		chrom1 = c1;
		start1 = s1;
		stop1 = st1;
		chrom2 = c2;
		start2 = s2;
		stop2 = st2;
		containsThree = false;
	}
	public void setChrom1(String c1){chrom1 = c1;}
	public String getChrom1(){return chrom1;}
	public void setStart1(int s1){start1 = s1;}
	public int getStart1(){return start1;}
	public void setStop1(int st1){stop1 = st1;}
	public int getStop1(){return stop1;}
	
	public void setChrom2(String c2){chrom2 = c2;}
	public String getChrom2(){return chrom2;}
	public void setStart2(int s2){start2 = s2;}
	public int getStart2(){return start2;}
	public void setStop2(int st2){stop2 = st2;}
	public int getStop2(){return stop2;}
	
	public void setChrom3(String c3){chrom3 = c3;}
	public String getChrom3(){return chrom3;}
	public void setStart3(int s3){start3 = s3;}
	public int getStart3(){return start3;}
	public void setStop3(int st3){stop3 = st3;}
	public int getStop3(){return stop3;}
	
	public boolean getContainsThree(){return containsThree;}
	
}
