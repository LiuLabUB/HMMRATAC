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
public class StateNode2 extends StateNode{

	
	private double emissionProb;
	private double transitionProb;
	private double zscore;
	public StateNode2(String c, int p, int s,double ep,double tp) {
		super(c, p, s);
		emissionProb = ep;
		transitionProb = tp;
	}
	public StateNode2(String c,int p,int s,double ep,double tp,double z){
		super(c, p, s);
		emissionProb = ep;
		transitionProb = tp;
		zscore = z;
	}
	
	public void setEP(double ep){
		emissionProb = ep;
	}
	public double getEP(){
		return emissionProb;
	}
	public void setTP(double tp){
		transitionProb = tp;
	}
	public double getTP(){
		return transitionProb;
	}
	public void setZScore(double z){zscore = z;}
	public double getZScore(){return zscore;}
}
