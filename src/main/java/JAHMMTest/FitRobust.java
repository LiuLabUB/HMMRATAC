package JAHMMTest;

import java.util.Arrays;

public class FitRobust {

	
	private double[][] covariance;
	
	
	public FitRobust(double[][] c){
		covariance = c;
		fit();
	}
	public double[][] getCovariance(){
		return covariance;
	}
	private void fit(){
		
		double[] m=new double [dimension(covariance)];
        for (int r = 0; r < dimension(covariance); r++){
                
                m[r]=covariance[r][r];
                        
                
        }
         Arrays.sort(m);
         double max = (m[m.length-1])/1000000000;
        
        
                
        
        while (!positiveDefined(covariance)){
                
                for (int r = 0; r < dimension(covariance); r++)
                        
                                covariance[r][r] = covariance[r][r]+max+0.00000000000000001;
                                        
        }
	}
	private boolean positiveDefined(double[][] m){
		double[][] l = new double[dimension(m)][dimension(m)];
		for (int j = 0; j < dimension(m); j++)
		{
			double[] lj = l[j];
			double d = 0.;
			
			for (int k = 0; k < j; k++) {
				double[] lk = l[k];
				double s = 0.;
				
				for (int i = 0; i < k; i++)
					s += lk[i] * lj[i];
				
				lj[k] = s = (m[j][k] - s) / l[k][k];
				d = d + s * s;
			}
			
			if ((d = m[j][j] - d) <= 0.)
				return false;
			
			l[j][j] = Math.sqrt(d);
			for (int k = j+1; k < dimension(m); k++)
				l[j][k] = 0.;
		}
		return true;
	}
	private int dimension(double[][] c){
		return c[0].length;
	}
}
