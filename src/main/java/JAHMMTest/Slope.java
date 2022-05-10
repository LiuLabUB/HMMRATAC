package JAHMMTest;

public class Slope {
	
	
	public Slope(double[] data, int w) {
		
	}
	
	public static double[] build(double[] data, int window) {
		double[] sloped = new double[data.length];
		int halfWidth = window / 2;
		for (int i = halfWidth; i < data.length - halfWidth; i++) {
			double start = data[i - halfWidth];
			double stop = data[i + halfWidth];
			double rise = stop - start;
			double run = window;
			if (window % 2 == 0) {
				run++;
			}
			sloped[i] = rise / run;
		}
		return sloped;
	}
}
