package stats;

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
