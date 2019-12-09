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
public class ZScore {
	
	private double[] _array;
	private double[] _transformedArray;
	private double _mean;
	
	public ZScore(double[] array){
		_array = array;
		_transformedArray = new double[_array.length];
		_mean = 0.0;
	}
	
	public void ZTransform(){
		double ave = 0.0;
		for (int i = 0; i < _array.length;i+=1){
			if (!Double.isNaN(_array[i])){
				if (_array[i] != 0.0){
					ave += _array[i];
				}
				else{
					ave+=1;//this applies only to Ratio program
				}
				
			}
			else{
				ave+=1;//note this else statement only applies to Ratio program.
			}
			
		}
		//the issue is that some of the values are Nan. when ave += Nan then ave = Nan
		
		ave /= (double)_array.length;
		_mean = ave;
		double var = 0.0;
		for (int i = 0; i < _array.length;i+=1){
			if (!Double.isNaN(_array[i])){
				if (_array[i] != 0.0){
					var += (_array[i]-ave)*(_array[i]-ave);
				}
				else{
					var += (1-ave)*(1-ave);
				}
			}
			else{
				var += (1-ave)*(1-ave);
			}
		}
		var /= (double)_array.length;
	
		double stdDev = Math.sqrt(var);
		//System.out.println("Mean"+"\t"+ave+"\t"+"Var"+"\t"+var);
		for (int i = 0; i < _array.length;i++){
			_transformedArray[i] = (_array[i] - ave) / stdDev;
		}
		
		
			
	}
	public double getMean(){
		return _mean;
	}
	
	public double[] getTransformedArray(){
		return _transformedArray;
	}

}
