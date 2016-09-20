package stats;

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
