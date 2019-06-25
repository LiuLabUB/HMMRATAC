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
public class GaussianKernel {
	
	private double[] _smoothed;
	private double[] _list;
	private int _degree;
	
	public GaussianKernel(double[] list,int degree){
		_degree = degree;
		_list = list;
		
		_smoothed = new double [list.length-(degree*2-1)];
	}
	
	public void smoothData(){
		double [] list = _list;
		int degree = _degree;
		
		int window = degree*2-1;
		double [] weight = new double[window];
		for (int i = 0;i < weight.length;i++){
			weight[i]++;
		}
		for (int i = 0;i < window;i++){
			float a = i-degree+1;
			float frac = a/(float) window;
			float first = 4*frac;
			float second = (float) Math.pow(first,2);
			float third = (float) Math.exp(second);
			double gauss = 1/third;
			weight[i] = gauss;
		}
		double sumweight = 0.0;
		for (int i = 0;i < weight.length;i++){
			sumweight += weight[i];
			
		}
		
		
		for (int i = 0;i < _smoothed.length;i++){
			
			double[] templist = new double[window];
			double sumtemplist=0.0;
			int counter=0;
			for (int a = i; a < i+window;a++){
				
				double wt = weight[counter];
				templist[counter]=list[a] * wt;
				sumtemplist+=templist[counter];
				counter++;
				
			}
			_smoothed[i]=sumtemplist/sumweight;
			
		}
		
		
	}
	
	public double[] getSmoothedData(){
		return _smoothed;
	}

}
