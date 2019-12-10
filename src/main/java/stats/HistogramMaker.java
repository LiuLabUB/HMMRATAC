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

public class HistogramMaker {
	private int[] _histogram;
	private double _small;
	private double _large;
	private double[] _peaks;
	private int _numberbins;
	
	public HistogramMaker(double small,double large,double[] Peaks,int numberbins){
		_small = small;
		_large = large;
		_peaks = Peaks;
		_numberbins = numberbins;
		
	}
	
	public void makeHistogram(){
		
		double range = _large - _small;
		double bins = range/_numberbins;
		_histogram = new int[_numberbins+1];
		
		int counter=0;
		for (double x = _small;x < _large;x+=bins){
			counter++;
			double maxinBin = x + bins;
			
			for (int i = 0;i < _peaks.length;i++){
				double score = _peaks[i];
				if (score >= x){
					if (score < maxinBin){
						_histogram[counter]++;
					}
				}
			}
		}	
		
	}
	
	public int[] getHistogram(){
		return _histogram;
	}

}
