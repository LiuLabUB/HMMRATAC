package stats;


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
