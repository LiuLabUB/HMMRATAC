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
public class DENode extends TagNode{

	private double average1;
	private double average2;
	private double pvalue;
	
	public DENode(String chr,int start, int stop,double ave1,double ave2,double p){
		super(chr,start,stop);
		average1 = ave1;
		average2 = ave2;
		pvalue = p;
	}
	public void setAve1(double ave1){
		average1 = ave1;
	}
	public double getAve1(){
		return average1;
	}
	public void setAve2(double ave2){
		average2 = ave2;
		
	}
	public double getAve2(){
		return average2;
	}
	public void setPValue(double p){
		pvalue = p;
	}
	public double getPValue(){
		return pvalue;
	}
	
}
