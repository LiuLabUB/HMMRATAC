package stats;
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
public class PearsonCorr {
	private double _corr = 0.0;
	
	public PearsonCorr(){

	}
	
	public PearsonCorr(double[] smoothed1,double[] smoothed2){
		int n = smoothed1.length +1;
		double totalone = 0.0;
		double totaltwo = 0.0;
		double totalonesquared = 0.0;
		double totaltwosquared = 0.0;
		double totalonetwo = 0.0;
		double totaloneone = 0.0;
		double totaltwotwo = 0.0;
		
		for (int a = 0; a < smoothed1.length;a++){
			double score1 = smoothed1[a];
			double score2 = smoothed2[a];
			totalone += score1;
			totaltwo += score2;
			totaloneone += score1*score1;
			totaltwotwo += score2*score2;
			totalonetwo += score1*score2;
		}
		totalonesquared = totalone*totalone;
		totaltwosquared = totaltwo*totaltwo;
		double topright = totalone*totaltwo;
		double topleft = n*totalonetwo;
		double top = topleft - topright;
		double bottomlefta = n*totaloneone;
		double bottomleftb = bottomlefta-totalonesquared;
		double bottomleft = Math.sqrt(bottomleftb);
		double bottomrighta = n*totaltwotwo;
		double bottomrightb = bottomrighta - totaltwosquared;
		double bottomright = Math.sqrt(bottomrightb);
		double bottom = bottomleft*bottomright;
		if (bottom > 0.0){
			_corr = top/bottom;
		}
		else {_corr = 0.0;}
		
	}
	

	public double getPearsonCorr(){
		return _corr;
	}

}
