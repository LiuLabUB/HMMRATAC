package JAHMMTest;

import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;


/**
 * This class can be used to compute the probability of a given observations
 * sequence for a given HMM.  Once the probability has been computed, the
 * object holds various information such as the <i>alpha</i> (and possibly
 * <i>beta</i>) array, as described in <i>Rabiner</i> and <i>Juang</i>.
 * <p>
 * Computing the <i>beta</i> array requires a O(1) access time to the 
 * observation sequence to get a theoretically optimal performance.
 */
public class ForwardBackwardCalculator_1
{       
        /**
         * Flags used to explain how the observation sequence probability
         * should be computed (either forward, using the alpha array, or backward,
         * using the beta array).
         */
        public static enum Computation { ALPHA, BETA };
        
        
        /* alpha[t][i] = P(O(1), O(2),..., O(t+1), i(t+1) = i+1 | hmm), that is the
         probability of the beginning of the state sequence (up to time t+1)
         with the (t+1)th state being i+1. */
        protected double[][] alpha = null;
        protected double[][] beta = null;
        protected double probability;
        
        
        protected ForwardBackwardCalculator_1()
        {
        };
        
        
        /**
         * Computes the probability of occurence of an observation sequence
         * given a Hidden Markov Model.
         *
         * @param hmm A Hidden Markov Model;
         * @param oseq An observation sequence.
         * @param flags How the computation should be done. See the
         *              {@link Computation Computation} enum.
         */
        public <O extends Observation>
        ForwardBackwardCalculator_1(List<? extends O> oseq,
                        Hmm<O> hmm, EnumSet<Computation> flags)
        {
        	System.out.println("Unscaled FB");
                if (oseq.isEmpty())
                        throw new IllegalArgumentException("Invalid empty sequence");
                
                if (flags.contains(Computation.ALPHA))
                        computeAlpha(hmm, oseq);
                
                if (flags.contains(Computation.BETA))
                        computeBeta(hmm, oseq);
                
                computeProbability(oseq, hmm, flags);
        }
        
        
        /**
         * Computes the probability of occurence of an observation sequence
         * given a Hidden Markov Model.  This computation computes the
         * <code>alpha</code> array as a side effect.
         * @see #ForwardBackwardCalculator(List, Hmm, EnumSet)
         */
        public <O extends Observation>
        ForwardBackwardCalculator_1(List<? extends O> oseq, Hmm<O> hmm)
        {
                this(oseq, hmm, EnumSet.of(Computation.ALPHA));
                
        }
        
        
        /* Computes the content of the alpha array */
        protected <O extends Observation> void
        computeAlpha(Hmm<? super O> hmm, List<O> oseq)
        {
                alpha = new double[oseq.size()][hmm.nbStates()];
                
                for (int i = 0; i < hmm.nbStates(); i++)
                        computeAlphaInit(hmm, oseq.get(0), i);
                
                Iterator<O> seqIterator = oseq.iterator();
                if (seqIterator.hasNext())
                        seqIterator.next();
                
                for (int t = 1; t < oseq.size(); t++) {
                        O observation = seqIterator.next();
                        
                        for (int i = 0; i < hmm.nbStates(); i++)
                                computeAlphaStep(hmm, observation, t, i);
                }
        }
        
        
        /* Computes alpha[0][i] */
        protected <O extends Observation> void
        computeAlphaInit(Hmm<? super O> hmm, O o, int i)
        {
                alpha[0][i] = hmm.getPi(i) * hmm.getOpdf(i).probability(o);
                if (alpha[0][i] == Double.NaN){System.out.println(alpha[0][i]);}
                //TEST OUT STATEMENT
                //System.out.println(alpha[0][i]);
                //System.out.println(hmm.getPi(i));
                //System.out.println(hmm.getOpdf(i).probability(o));
        }
        
        
        /* Computes alpha[t][j] (t > 0) */
        protected <O extends Observation> void 
        computeAlphaStep(Hmm<? super O> hmm, O o, int t, int j)
        {
                double sum = 0.;
                
                for (int i = 0; i < hmm.nbStates(); i++){
                        sum += alpha[t-1][i] * hmm.getAij(i, j);
                        //TEST OUT STATEMENT
                        //System.out.println(hmm.getAij(i, j));
                        //System.out.println(alpha[t-1][i]);
                }

                alpha[t][j] = sum * hmm.getOpdf(j).probability(o);
                //TEST OUT STATEMENT
              // if (Double.isNaN(alpha[t][j])){System.out.println("Alpha"+"\t"+alpha[t][j]+"\tprobability\t"+hmm.getOpdf(j).probability(o)
            	//	   +"\tstate\t"+j+"\tobservation\t"+o.toString()+"\t"+o);}
                //System.out.println(hmm.getOpdf(j).probability(o));
                //System.out.println(sum);
        }
        
        
        /* Computes the content of the beta array.  Needs a O(1) access time
         to the elements of oseq to get a theoretically optimal algorithm. */
        protected <O extends Observation> void 
        computeBeta(Hmm<? super O> hmm, List<O> oseq)
        {
                beta = new double[oseq.size()][hmm.nbStates()];
                
                for (int i = 0; i < hmm.nbStates(); i++)
                        beta[oseq.size()-1][i] = 1.;
                
                for (int t = oseq.size()-2; t >= 0; t--)
                        for (int i = 0; i < hmm.nbStates(); i++)
                                computeBetaStep(hmm, oseq.get(t+1), t, i);
        }
        
        
        /* Computes beta[t][i] (t < obs. seq.le length - 1) */
        protected <O extends Observation> void 
        computeBetaStep(Hmm<? super O> hmm, O o, int t, int i)
        {
                double sum = 0.;
                
                for (int j = 0; j < hmm.nbStates(); j++)
                        sum += beta[t+1][j] * hmm.getAij(i, j) * 
                        hmm.getOpdf(j).probability(o);
                
                beta[t][i] = sum;
                if (beta[t][i] == Double.NaN){System.out.println(beta[t][i]);}
        }
        
        
        /**
         * Returns an element of the <i>alpha</i> array.
         * 
         * @param t The temporal argument of the array (positive but strictly
         *          smaller than the length of the sequence that helped generating
         *          the array).
         * @param i A state index of the HMM that helped generating the array.
         * @throws {@link UnsupportedOperationException 
         *          UnsupportedOperationException} if alpha array has not been
         *          computed.
         * @return The <i>alpha</i> array (t, i) element.
         */ 
        public double alphaElement(int t, int i)
        {
                if (alpha == null)
                        throw new UnsupportedOperationException("Alpha array has not " +
                                        "been computed");
                
                return alpha[t][i];
        }
        
        
        /**
         * Returns an element of the <i>beta</i> array.
         * 
         * @param t The temporal argument of the array (positive but smaller than
         *          the length of the sequence that helped generating the array).
         * @param i A state index of the HMM that helped generating the array.
         * @throws {@link UnsupportedOperationException 
         *          UnsupportedOperationException} if beta array has not been
         *          computed.
         * @return The <i>beta</i> beta (t, i) element.
         */ 
        public double betaElement(int t, int i)
        {
                if (beta == null)
                        throw new UnsupportedOperationException("Beta array has not " +
                                        "been computed");
                
                return beta[t][i];
        }
        
        
        private <O extends Observation> void 
        computeProbability(List<O> oseq, Hmm<? super O> hmm, 
                        EnumSet<Computation> flags)
        {
                probability = 0.;
                
                if (flags.contains(Computation.ALPHA))
                        for (int i = 0; i < hmm.nbStates(); i++){ 
                                probability += alpha[oseq.size()-1][i];
                               // System.out.println(probability);
                        }
                else
                        for (int i = 0; i < hmm.nbStates(); i++)
                                probability += 
                                        hmm.getPi(i) * 
                                        hmm.getOpdf(i).probability(oseq.get(0)) * beta[0][i];
        }
        
        
        /**
         * Return the probability of the sequence that generated this object.
         * For long sequences, this probability might be very small, or even
         * meaningless because of underflows.
         *
         * @return The probability of the sequence of interest.
         */
        public double probability()
        {
        	//TEST OUT STATEMENT
        	//System.out.println("Program is using correct ForwardBackward Class!!");
                return probability;
        }
}