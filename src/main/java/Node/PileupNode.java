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
public class PileupNode {
	
	
	private int _base;
	private double _score;
	private int _score2;
	/**
	 * Constructor for creating new PileupNode
	 * @param base an integer representing the base position
	 * @param score a double representing the base's score
	 */
	public PileupNode(int base, double score){
	//	_chrom = chr;
		_base = base;
		_score = score;
	}
	/**
	 * Default Constructor
	 */
	public PileupNode(){
		
	}
	/**
	 * Constructor for creating new PileupNode
	 * @param base an integer representing the base position
	 * @param score a integer representing the base's score
	 */
	public PileupNode(int base,int score){
		_base = base;
		_score2 = score;
	}
	/**
	 * Access the base position
	 * @return an integer representing the base position
	 */
	public int getBase(){
		return _base;
	}
	/**
	 * Set the base position
	 * @param base an integer representing the base position
	 */
	public void setBase(int base){
		_base = base;
	}
	/**
	 * Access the score
	 * @return a double representing the score at the position
	 */
	public double getScore(){
		return _score;
	}
	/**
	 * Set the score
	 * @param score a double representing the positions score
	 */
	public void setScore(double score){
		_score = score;
	}
	/**
	 * Access the score
	 * @return an integer representing the score at the position
	 */
	public int getScore2(){
		return _score2;
	}
	/**
	 * Set the score
	 * @param score an integer representing the positions score
	 */
	public void setScore2(int score){
		_score2 = score;
	}
}
