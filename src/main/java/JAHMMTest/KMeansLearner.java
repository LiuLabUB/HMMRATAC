package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the K-Means learning algorithm.
 */
public class KMeansLearner<O extends Observation & CentroidFactory<? super O>>
{       
        private Clusters<O> clusters;
        private int nbStates;
        private List<? extends List<? extends O>> obsSeqs;
        private OpdfFactory<? extends Opdf<O>> opdfFactory;
        private boolean terminated;
        
        
        /**
         * Initializes a K-Means algorithm implementation.  This algorithm
         * finds a HMM that models a set of observation sequences.
         *
         * @param nbStates  The number of states the resulting HMM will be made of.
         * @param opdfFactory A class that builds the observation probability
         *                    distributions associated to the states of the HMM.
         * @param sequences A vector of observation sequences.  Each observation
         *                sequences is a vector of
         *                {@link be.ac.ulg.montefiore.run.jahmm.Observation
         *                observations} compatible with the 
         *                {@link be.ac.ulg.montefiore.run.jahmm.CentroidFactory
         *                k-means algorithm}.
         */
        public KMeansLearner(int nbStates,
                        OpdfFactory<? extends Opdf<O>> opdfFactory,
                        List<? extends List<? extends O>> sequences)
        {       
                this.obsSeqs = sequences;
                this.opdfFactory = opdfFactory;
                this.nbStates = nbStates;
                
                List<? extends O> observations = flat(sequences);
                clusters = new Clusters<O>(nbStates, observations);
                terminated = false;
        }
        
        
        /**
         * Performs one iteration of the K-Means algorithm.
         * In one iteration, a new HMM is computed using the current clusters, and
         * the clusters are re-estimated using this HMM.
         *
         * @return A new, updated HMM.
         */
        public Hmm<O> iterate()
        {       
                Hmm<O> hmm = new Hmm<O>(nbStates, opdfFactory);
                
                learnPi(hmm);
                learnAij(hmm);
                learnOpdf(hmm);
                
                terminated = optimizeCluster(hmm);
                
                return hmm;
        }
        
        
        /**
         * Returns <code>true</code> if the algorithm has reached a fix point,
         * else returns <code>false</code>.
         */
        public boolean isTerminated()
        {
                return terminated;
        }
        
        
        /**
         * Does iterations of the K-Means algorithm until a fix point is reached.
         * 
         * @return The HMM that best matches the set of observation sequences given
         *         (according to the K-Means algorithm).
         */
        public Hmm<O> learn()
        {       
                Hmm<O> hmm;
                
                do 
                        hmm = iterate();
                while(!isTerminated());
                
                return hmm;
        }
        
        
        private void learnPi(Hmm<?> hmm)
        {       
                double[] pi = new double[nbStates];
                
                for (int i = 0; i < nbStates; i++)
                        pi[i] = 0.;
                
                for (List<? extends O> sequence : obsSeqs)
                        pi[clusters.clusterNb(sequence.get(0))]++;
                
                for (int i = 0; i < nbStates; i++)
                        hmm.setPi(i, pi[i] / obsSeqs.size());
        }
        
        
        private void learnAij(Hmm<O> hmm)
        {       
                for (int i = 0; i < hmm.nbStates(); i++)
                        for (int j = 0; j < hmm.nbStates(); j++)
                                hmm.setAij(i, j, 0.);
                
                for (List<? extends O> obsSeq : obsSeqs) {
                        if (obsSeq.size() < 2)
                                continue;
                        
                        int first_state;
                        int second_state = clusters.clusterNb(obsSeq.get(0));
                        for (int i = 1; i < obsSeq.size(); i++) {
                                first_state = second_state;
                                second_state =
                                        clusters.clusterNb(obsSeq.get(i));
                                
                                hmm.setAij(first_state, second_state,
                                                hmm.getAij(first_state, second_state)+1.);
                        }
                }
                
                /* Normalize Aij array */
                for (int i = 0; i < hmm.nbStates(); i++) {
                        double sum = 0;
                        
                        for (int j = 0; j < hmm.nbStates(); j++)
                                sum += hmm.getAij(i, j);
                        
                        if (sum == 0.)
                                for (int j = 0; j < hmm.nbStates(); j++) 
                                        hmm.setAij(i, j, 1. / hmm.nbStates());     // Arbitrarily
                        else
                                for (int j = 0; j < hmm.nbStates(); j++)
                                        hmm.setAij(i, j, hmm.getAij(i, j) / sum);
                }
        }
        
        
        private void learnOpdf(Hmm<O> hmm)
        {
                for (int i = 0; i < hmm.nbStates(); i++) {
                        Collection<O> clusterObservations = clusters.cluster(i);
                        
                        if (clusterObservations.isEmpty())
                                hmm.setOpdf(i, opdfFactory.factor());
                        else
                                hmm.getOpdf(i).fit(clusterObservations);
                }
        }
        
        
        /* Return true if no modification */
        private boolean optimizeCluster(Hmm<O> hmm)
        {       
                boolean modif = false;
                
                for (List<? extends O> obsSeq : obsSeqs) {
                        ViterbiCalculator vc = new ViterbiCalculator(obsSeq, hmm);
                        int states[] = vc.stateSequence();
                        
                        for (int i = 0; i < states.length; i++) {
                                O o = obsSeq.get(i);
                                
                                if (clusters.clusterNb(o) != states[i]) {
                                        modif = true;
                                        clusters.remove(o, clusters.clusterNb(o));
                                        clusters.put(o, states[i]);
                                }
                        }
                }
                
                return !modif;
        }
        
        
        public static <T> List<T> flat(List<? extends List<? extends T>> lists)
        {       
                List<T> v = new ArrayList<T>();
                
                for (List<? extends T> list : lists)
                        v.addAll(list);
                
                return v;
        }
}


/*
 * This class holds the matching between observations and clusters.
 */
class Clusters<O extends CentroidFactory<? super O>>
{       
        class Value
        {
                private int clusterNb;
                
                Value(int clusterNb)
                {
                        this.clusterNb = clusterNb;
                }
                
                void setClusterNb(int clusterNb)
                {
                        this.clusterNb = clusterNb;
                }
                
                int getClusterNb()
                {
                        return clusterNb;
                }
        }
        
        
        private Hashtable<O,Value> clustersHash;
        private ArrayList<Collection<O>> clusters;
        
        
        public Clusters(int k, List<? extends O> observations)
        {
                
                clustersHash = new Hashtable<O,Value>();
                clusters = new ArrayList<Collection<O>>();
                
                KMeansCalculator<O> kmc = new KMeansCalculator<O>(k, observations);
                
                for (int i = 0; i < k; i++) {
                        Collection<O> cluster = kmc.cluster(i);
                        clusters.add(cluster);
                        
                        for (O element : cluster) 
                                clustersHash.put(element, new Value(i));
                }
        }
        
        
        public boolean isInCluster(Observation o, int clusterNb)
        {
                return clusterNb(o) == clusterNb;
        }
        
        
        public int clusterNb(Observation o)
        {
                return clustersHash.get(o).getClusterNb();
        }
        
        
        public Collection<O> cluster(int clusterNb)
        {
                return clusters.get(clusterNb);
        }
        
        
        public void remove(Observation o, int clusterNb)
        {
                clustersHash.get(o).setClusterNb(-1);
                clusters.get(clusterNb).remove(o);
        }
        
        
        public void put(O o, int clusterNb)
        {
                clustersHash.get(o).setClusterNb(clusterNb);
                clusters.get(clusterNb).add(o);
        }
}/*
Hide details
Change log
r34 by jm.francois on Mar 1, 2009   Diff
removed license, renewed build
Go to: 	
Older revisions
 r19 by jm.francois on Feb 15, 2009   Diff 
All revisions of this file
File info
Size: 6193 bytes, 278 lines
View raw file
File properties
copyright
Copyright (c) 2004-2009, Jean-Marc Francois*/
