package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the Baum-Welch learning algorithm.  This algorithm
 * finds a HMM that models a set of observation sequences.
 */
public class BaumWelchLearner
{       
        /**
         * Number of iterations performed by the {@link #learn} method.
         */
        private int nbIterations = 9;
        @SuppressWarnings("unused")
		private int constant = 1;
        
        /**
         * Initializes a Baum-Welch instance. 
         */
        protected ArrayList<Double> pts = new ArrayList<Double>();
        public BaumWelchLearner()
        {
        }
        public ArrayList<Double> getPoints(){return pts;}
        
        /**
         * Performs one iteration of the Baum-Welch algorithm.
         * In one iteration, a new HMM is computed using a previously estimated
         * HMM.
         *
         * @param hmm A previously estimated HMM.
         * @param sequences The observation sequences on which the learning is
         *         based.  Each sequence must have a length higher or equal to
         *         2.
         * @return A new, updated HMM.
         */
        public <O extends Observation> Hmm<O>
        iterate(Hmm<O> hmm, List<? extends List<? extends O>> sequences)
        {               
                Hmm<O> nhmm;
                try {
                        nhmm = hmm.clone();
                } catch(CloneNotSupportedException e) {
                        throw new InternalError();
                }
                        
                /* gamma and xi arrays are those defined by Rabiner and Juang */
                /* allGamma[n] = gamma array associated to observation sequence n */
                double allGamma[][][] = new double[sequences.size()][][];
                
                /* a[i][j] = aijNum[i][j] / aijDen[i]
                 * aijDen[i] = expected number of transitions from state i
                 * aijNum[i][j] = expected number of transitions from state i to j
                 */
                double aijNum[][] = new double[hmm.nbStates()][hmm.nbStates()];
                double aijDen[] = new double[hmm.nbStates()];
                
                Arrays.fill(aijDen, 0.);
                for (int i = 0; i < hmm.nbStates(); i++)
                        Arrays.fill(aijNum[i], 0.);
                
                int g = 0;
                for (List<? extends O> obsSeq : sequences) {        
                        ForwardBackwardCalculator_1 fbc = 
                                generateForwardBackwardCalculator(obsSeq, hmm);
                        
                        double xi[][][] = estimateXi(obsSeq, fbc, hmm);
                        double gamma[][] = allGamma[g++] = estimateGamma(xi, fbc);
                        
                        for (int i = 0; i < hmm.nbStates(); i++)
                                for (int t = 0; t < obsSeq.size() - 1; t++) {
                                        aijDen[i] += gamma[t][i];
                                        
                                        for (int j = 0; j < hmm.nbStates(); j++){
                                                aijNum[i][j] += xi[t][i][j];
                                                //TEST OUT STATEMENT
                                               // System.out.println(aijNum[i][j]);
                                        }
                                }
                }
                
                for (int i = 0; i < hmm.nbStates(); i++) {
                        if (aijDen[i] == 0.) // State i is not reachable
                                for (int j = 0; j < hmm.nbStates(); j++)
                                        nhmm.setAij(i, j, hmm.getAij(i, j));
                        else
                                for (int j = 0; j < hmm.nbStates(); j++)
                                        nhmm.setAij(i, j, aijNum[i][j] / aijDen[i]);
                }
                
                /* pi computation */
                for (int i = 0; i < hmm.nbStates(); i++)
                        nhmm.setPi(i, 0.);
                
                for (int o = 0; o < sequences.size(); o++)
                        for (int i = 0; i < hmm.nbStates(); i++){
                                nhmm.setPi(i,nhmm.getPi(i) + allGamma[o][0][i] / sequences.size());
                                //nhmm.setPi(i, nhmm.getPi(i));
                                //TEST OUT STATEMENT
                                //System.out.println(allGamma[o][0][i]);
                        }
                
                
                /* pdfs computation */
                //pts = new ArrayList<Double>();
                for (int i = 0; i < hmm.nbStates(); i++) {
                        List<O> observations = flat(sequences);
                        double[] weights = new double[observations.size()];
                        double sum = 0.;
                        int j = 0;
                        
                        int o = 0;
                        //Test out statement
                       // System.out.println(sequences.size());
                        for (List<? extends O> obsSeq : sequences) {
                        	//TEST OUT STATEMENT
                        	//System.out.println(obsSeq.size());
                                for (int t = 0; t < obsSeq.size(); t++, j++){
                                        sum += weights[j] = allGamma[o][t][i];
                                        //TEST OUT STATEMENT
                                       //System.out.println(allGamma[o][t][i]);
                                }
                                o++;
                        }
                        
                        for (j--; j >= 0; j--){
                                weights[j] /= sum;
                               // pts.add(weights[j]);
                                //System.out.println(weights[j]);
                        }
                        Opdf<O> opdf = nhmm.getOpdf(i);
                       
                        //added lines that enforce a diagonal cov matrix
                        /*
                        OpdfMultiGaussian temp = (OpdfMultiGaussian) opdf;
                        OpdfMultiGaussian_2 pdf = new OpdfMultiGaussian_2(temp.mean(),temp.covariance());
                        pdf.fit((List<ObservationVector>)observations, weights);
                        temp = new OpdfMultiGaussian(pdf.mean(),pdf.covariance());
                        nhmm.setOpdf(i, (Opdf<O>) temp);
                        */
                        
                        //Below line allows for non-diagonal cov matrix
                        opdf.fit(observations,weights);
                        
                }
                // Last step out statement:
                //System.out.println(nhmm.toString());
                return nhmm;
        }
        public void setConstant(int c){constant = c;}
        
        protected <O extends Observation> ForwardBackwardCalculator_1
        generateForwardBackwardCalculator(List<? extends O> sequence, Hmm<O> hmm)
        {       
                return new ForwardBackwardCalculator_1(sequence, hmm, 
                                EnumSet.allOf(ForwardBackwardCalculator_1.Computation.class));
        }
        
        
        /**
         * Does a fixed number of iterations (see {@link #getNbIterations}) of the
         * Baum-Welch algorithm.
         * 
         * @param initialHmm An initial estimation of the expected HMM.  This
         *         estimate is critical as the Baum-Welch algorithm only find
         *         local minima of its likelihood function.
         * @param sequences The observation sequences on which the learning is
         *         based.  Each sequence must have a length higher or equal to 2.
         * @return The HMM that best matches the set of observation sequences given
         *         (according to the Baum-Welch algorithm).
         */
        public <O extends Observation> Hmm<O>
        learn(Hmm<O> initialHmm, List<? extends List<? extends O>> sequences)
        {
                Hmm<O> hmm = initialHmm;
                
                for (int i = 0; i < nbIterations; i++)
                        hmm = iterate(hmm, sequences);
                
                return hmm;
        }
        
        
        protected <O extends Observation> double[][][]
        estimateXi(List<? extends O> sequence, ForwardBackwardCalculator_1 fbc,
                        Hmm<O> hmm)
        {       
                if (sequence.size() <= 1)
                        throw new IllegalArgumentException("Observation sequence too " +
                        "short");
                
                double xi[][][] = 
                        new double[sequence.size()-1][hmm.nbStates()][hmm.nbStates()];
                double probability = fbc.probability();
              //  System.out.println("Xi Method, unscaled");
                //TEST OUT STATEMENT
                //System.out.println(probability);
                Iterator<? extends O> seqIterator = sequence.iterator();
                seqIterator.next();
                
                //added line
               // pts = new double[sequence.size()-1][hmm.nbStates()];
                
                
                for (int t = 0; t < sequence.size() - 1; t++) {
                        O o = seqIterator.next();
                       // System.out.println("First Loop Xi unscaled");
                        for (int i = 0; i < hmm.nbStates(); i++)
                                for (int j = 0; j < hmm.nbStates(); j++){
                                		//System.out.println("Second Loop Xi unscaled");
                                        xi[t][i][j] = fbc.alphaElement(t, i) *
                                        hmm.getAij(i, j) *
                                        hmm.getOpdf(j).probability(o) *
                                        fbc.betaElement(t+1, j) / probability;
                                        
                                        double value = fbc.alphaElement(t, i) *
                                        		hmm.getAij(i, j) *
                                        		hmm.getOpdf(j).probability(o) *
                                        		fbc.betaElement(t+1, j);
                                       if (Double.isNaN(value)){
                                    	   /*
                                    	   System.out.println("xi[t][i][j]"+"\t"+value);
                                    	   System.out.println("aplha"+"\t"+fbc.alphaElement(t, i));
                                    	   System.out.println("beta"+"\t"+fbc.betaElement(t+1, j));
                                    	   System.out.println("transition"+"\t"+hmm.getAij(i, j));
                                    	   System.out.println("emission\t"+hmm.getOpdf(j).probability(o));
                                    	   */
                                       }
                                        xi[t][i][j] = value;
                                       // pts[t][i] = fbc.alphaElement(t, i)
                                        		
                                        		;
                                        //Testing the probability
                                        /*OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(j);
                                        double[] means = pdf.mean();
                                        double[][] cov = pdf.covariance();
                                        MultivariateNormalDistribution normal = new MultivariateNormalDistribution(means,cov);
                                        ObservationVector obs = (ObservationVector) o;
                                        double[] val = obs.values();
                                        System.out.println("Apache density:"+"\t"+normal.density(val));
                                        */
                                        
                                        /*
                                         * Convert the multivariate dist into multiple monovariate gaussians,
                                         * calculate the probability of each function and multiply together
                                         */
                                        /*
                                        OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(j);
                                        double[] means = pdf.mean();
                                        double[][] cov = pdf.covariance();
                                        ObservationVector obs = (ObservationVector) o;
                                        double[] val = obs.values();
                                        double density = 1.0;
                                        for (int a = 0;a < means.length;a++){
                                        	double mean = means[a];
                                        	double var = cov[a][a];
                                        	double sd = Math.sqrt(var);
                                        	NormalDistribution normal = new NormalDistribution(mean,sd);
                                        	density *= normal.probability(val[a]-0.000000001,val[a]+0.000000001);
                                        }*/
                                        //if(density >= 1.0){
                                        	//System.out.println("Mixture Mono Gaussian Density:"+"\t"+density);
                                        //}
                                        //TEST OUT STATEMENT
                                        //System.out.println(hmm.getAij(i, j));
                                        //if (hmm.getOpdf(j).probability(o) >= 1.0){
                                        	//System.out.println("OPDF density"+"\t"+hmm.getOpdf(j).probability(o));
                                        //}
                                        //System.out.println(o.toString());
                                        //System.out.println("Sequence:"+"\t"+t+"\t"+"State:"+"\t"+i);
                                        //System.out.println(fbc.alphaElement(t, i));
                                        //System.out.println(fbc.betaElement(t+1, j));
                                        //System.out.println(probability);
                                }
                }
                
                return xi;
        }
        
        
        /* gamma[][] could be computed directly using the alpha and beta
         * arrays, but this (slower) method is prefered because it doesn't
         * change if the xi array has been scaled (and should be changed with
         * the scaled alpha and beta arrays).
         */
        protected double[][]
        estimateGamma(double[][][] xi, ForwardBackwardCalculator_1 fbc)
        {
                double[][] gamma = new double[xi.length + 1][xi[0].length];
                
                for (int t = 0; t < xi.length + 1; t++)
                        Arrays.fill(gamma[t], 0.);
                
                for (int t = 0; t < xi.length; t++)
                        for (int i = 0; i < xi[0].length; i++)
                                for (int j = 0; j < xi[0].length; j++){
                                        gamma[t][i] += xi[t][i][j];
                                        //TEST OUT STATEMENT
                                        /*
                                       if (Double.isNaN(gamma[t][i])){ System.out.println(gamma[t][i]);
                                       System.out.println(xi[t][i][j]);
                                       System.out.println("Xi"+t+"_"+i+"_"+j);}*/
                                }
                
                for (int j = 0; j < xi[0].length; j++)
                        for (int i = 0; i < xi[0].length; i++)
                                gamma[xi.length][j] += xi[xi.length - 1][i][j];
                
                return gamma;
        }
        
        
        /**
         * Returns the number of iterations performed by the {@link #learn} method.
         * 
         * @return The number of iterations performed.
         */
        public int getNbIterations()
        {
                return nbIterations;
        }
        
        
        /**
         * Sets the number of iterations performed by the {@link #learn} method.
         * 
         * @param nb The (positive) number of iterations to perform.
         */
        public void setNbIterations(int nb)
        {
                if (nb < 0)
                        throw new IllegalArgumentException("Positive number expected");
                
                nbIterations = nb;
        }
        public static <T> List<T> flat(List<? extends List<? extends T>> lists)
        {       
                List<T> v = new ArrayList<T>();
                
                for (List<? extends T> list : lists)
                        v.addAll(list);
                
                return v;
        }
}
