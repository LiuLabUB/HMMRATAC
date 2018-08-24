package JAHMMTest;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;

import be.ac.ulg.montefiore.run.distributions.MultiGaussianDistribution;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.Opdf;



/**
 * This class represents a multivariate gaussian distribution function.
 */
public class OpdfMultiGaussian_2
implements Opdf<ObservationVector>
{       
        private MultiGaussianDistribution distribution;
        
        
        /**
         * Builds a new gaussian probability distribution with zero mean and
         * identity covariance matrix.
         *
         * @param dimension The dimension of the vectors.
         */
        public OpdfMultiGaussian_2(int dimension)
        {
                distribution = new MultiGaussianDistribution(dimension);
        }
        
        
        /**
         * Builds a new gaussian probability distribution with a given mean and
         * covariance matrix.
         *
         * @param mean The distribution's mean.
         * @param covariance The distribution's covariance matrix.
         */
        public OpdfMultiGaussian_2(double[] mean, double[][] covariance)
        {               
                if (covariance.length == 0 || mean.length != covariance.length ||
                                covariance.length != covariance[0].length)
                        throw new IllegalArgumentException();
                
                distribution = new MultiGaussianDistribution(mean, covariance);
        }
        
        
        /**
         * Returns (a copy of) this distribution's mean vector.
         *
         * @return The mean vector.
         */
        public double[] mean()
        {
                return distribution.mean();
        }
        
        
        /**
         * Returns (a copy of) this distribution's covariance matrix.
         *
         * @return The covariance matrix.
         */
        public double[][] covariance()
        {
                return distribution.covariance();
        }
        
        
        /**
         * Returns the dimension of the vectors handled by this distribution.
         *
         * @return The dimension of the vectors handled by this distribution.
         */
        public int dimension()
        {
                return distribution.dimension();
        }
        
        
        public double probability(ObservationVector o)
        {
                if (o.dimension() != distribution.dimension())
                        throw new IllegalArgumentException("Vector has a wrong " +
                        "dimension");
                
                return distribution.probability(o.values());
        }
        
        
        public ObservationVector generate()
        {
                return new ObservationVector(distribution.generate());
        }
        
        
        public void fit(ObservationVector... oa)
        {
                fit(Arrays.asList(oa));
        }
        
        
        public void fit(Collection<? extends ObservationVector> co)
        {
                if (co.isEmpty())
                        throw new IllegalArgumentException("Empty observation set");
                
                double[] weights = new double[co.size()];
                Arrays.fill(weights, 1. / co.size());
                
                fit(co, weights);
        }
        
        
        public void fit(ObservationVector[] o, double[] weights)
        {
                fit(Arrays.asList(o), weights);
        }
        
        
        public void fit(Collection<? extends ObservationVector> co, 
                        double[] weights)
        {
                if (co.isEmpty() || co.size() != weights.length)
                        throw new IllegalArgumentException();
                
                // Compute mean
                double[] mean = new double[dimension()];
                for (int r = 0; r < dimension(); r++) {
                        int i = 0;
                        
                        for (ObservationVector o : co)
                                mean[r] += o.values()[r] * weights[i++];
                }
                
                // Compute covariance
                double[][] covariance = new double[dimension()][dimension()];
                int i = 0;
                for (ObservationVector o : co) {
                        double[] obs = o.values();
                        double[] omm = new double[obs.length];
                        
                        for (int j = 0; j < obs.length; j++)
                                omm[j] = obs[j] - mean[j];
                        
                        for (int r = 0; r < dimension(); r++)
                                for (int c = 0; c < dimension(); c++)
                                	if (r == c)//added line
                                        covariance[r][c] += omm[r] * omm[c] * weights[i];
                        
                        i++;
                }
                
                distribution = new MultiGaussianDistribution(mean, covariance);
        }
        
        
        public OpdfMultiGaussian_2 clone()
        {
                try {
                        return (OpdfMultiGaussian_2) super.clone();
                } catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
        }
        
        
        public String toString()
        {
                return toString(NumberFormat.getInstance());
        }
        
        
        public String toString(NumberFormat numberFormat)
        {
                String s = "Multi-variate Gaussian distribution --- Mean: [ ";
                double[] mean = distribution.mean();
                
                for (int i = 0; i < mean.length; i++)
                        s += numberFormat.format(mean[i]) + " ";
                
                return s + "]";
        }


        private static final long serialVersionUID = 1L;
}