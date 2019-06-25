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


//should re-do this class.  should make one basic class (for 6 arg constructor) and then make other classes that implement this one
//for other functions - for instance for node with cluster

public class ATACMatrixNode implements Comparable{

	
	private String _chrom;
	private int _pos;
	private double _enrich1;
	private double _enrich2;
	private double _enrich3;
	private double _enrich4;
	
	private double _windowSum1;
	private double _windowSum2;
	private double _windowSum3;
	
	private double _lamda1;
	private double _lamda2;
	private double _lamda3;
	
	private int _index;
	
	public ATACMatrixNode(String chrom,int pos, double enrich1,double enrich2,double enrich3,int index){
		_chrom = chrom;
		_pos = pos;
		_enrich1 = enrich1;
		_enrich2 = enrich2;
		_enrich3 = enrich3;
		_index = index;
	}
	
	public ATACMatrixNode(String chrom,int pos, double enrich1,double enrich2,double enrich3,double enrich4,int index){
		_chrom = chrom;
		_pos = pos;
		_enrich1 = enrich1;
		_enrich2 = enrich2;
		_enrich3 = enrich3;
		_enrich4 = enrich4;
		_index = index;
	}
	
	public ATACMatrixNode(String chrom,int pos,double enrich1,double enrich2,double enrich3,double winSum1,double winSum2,double winSum3,double lamda1,double lamda2,double lamda3){
		_chrom = chrom;
		_pos = pos;
		_enrich1 = enrich1;
		_enrich2 = enrich2;
		_enrich3 = enrich3;
		
		_windowSum1 = winSum1;
		_windowSum2 = winSum2;
		_windowSum3 = winSum3;
		
		_lamda1 = lamda1;
		_lamda2 = lamda2;
		_lamda3 = lamda3;
	}
	
	public void setIndex(int index){
		_index=index;
	}
	public int getIndex(){
		return _index;
	}
	
	public void setChrom(String chr){
		_chrom = chr;
	}
	public String getChrom(){
		return _chrom;
	}
	public void setPos(int pos){
		_pos = pos;
	}
	public int getPos(){
		return _pos;
	}
	public void setEnrich1(double enrich1){
		_enrich1 = enrich1;
	}
	public double getEnrich1(){
		return _enrich1;
	}
	public void setEnrich2(double enrich2){
		_enrich2 = enrich2;
	}
	public double getEnrich2(){
		return _enrich2;
	}
	public void setEnrich3(double enrich3){
		_enrich3 = enrich3;
	}
	public double getEnrich3(){
		return _enrich3;
	}
	public void setEnrich4(double enrich4){
		_enrich4 = enrich4;
	}
	public double getEnrich4(){
		return _enrich4;
	}
	
	public void setSum1(double winSum1){
		_windowSum1 = winSum1;
	}
	public double getSum1(){
		return _windowSum1;
	}
	public void setSum2(double winSum2){
		_windowSum2 = winSum2;
	}
	public double getSum2(){
		return _windowSum2;
	}
	public void setSum3(double winSum3){_windowSum3=winSum3;}
	public double getSum3(){return _windowSum3;}
	
	public void setLamda1(double lamda1){
		_lamda1 = lamda1;
	}
	public double getLamda1(){
		return _lamda1;
	}
	public void setLamda2(double lamda2){
		_lamda2 = lamda2;
	}
	public double getLamda2(){
		return _lamda2;
	}
	public void setLamda3(double local3){_lamda3=local3;}
	public double getLamda3(){return _lamda3;}
	
	public static Comparator<ATACMatrixNode> positionComparator = new Comparator<ATACMatrixNode>(){
		public int compare(ATACMatrixNode node1, ATACMatrixNode node2){
			double pos1 = node1.getPos();
			double pos2 = node2.getPos();
			if (pos1 < pos2) return -1;
			else return 1;
		
		}
	};

	

	@Override
	public int compareTo(Object o) {
		double pos1 = this.getPos();
		double pos2 = ((ATACMatrixNode) o).getPos();
		if (pos1 < pos2)return -1;
		else if (pos1 > pos2)return 1;
		else return 0;
	}

}
