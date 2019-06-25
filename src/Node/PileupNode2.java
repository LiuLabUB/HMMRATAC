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
public class PileupNode2 extends PileupNode{

	
	private String chrom;
	
	/**
	 * Constructor for creating PileupNode2 object
	 * @param base an integer representing the base position
	 * @param score an integer representing the position's score
	 * @param chr a String representing the chromosome name
	 */
	public PileupNode2(int base,int score,String chr){
		super(base,score);
		chrom = chr;
	}
	/**
	 * Constructor for creating PileupNode2 object
	 * @param base an integer representing the base position
	 * @param score a double representing the position's score
	 * @param chr a String representing the chromosome name
	 */
	public PileupNode2(int base,double score,String chr){
		super(base,score);
		chrom = chr;
	}
	
	
	/**
	 * Set the chromsome name
	 * @param chr a String representing the chromosome name
	 */
	public void setChrom(String chr){
		chrom = chr;
	}
	/**
	 * Access the chromosome name
	 * @return a String representing the chromosome na,e
	 */
	public String getChrom(){
		return chrom;
	}
	
}
