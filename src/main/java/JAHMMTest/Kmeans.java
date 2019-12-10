package JAHMMTest;


import java.util.ArrayList;
import java.util.Random;

import net.sf.javaml.clustering.Clusterer;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.distance.DistanceMeasure;
import net.sf.javaml.distance.EuclideanDistance;
import net.sf.javaml.tools.DatasetTools;

/**
 * Implements the K-means algorithms as described by Mac Queen in 1967.
 * 
 * <bibtex> J. B. MacQueen (1967): "Some Methods for classification and Analysis
 * of Multivariate Observations, Proceedings of 5-th Berkeley Symposium on
 * Mathematical Statistics and Probability", Berkeley, University of California
 * Press, 1:281-297 </bibtex>
 * 
 * 
 * @author Thomas Abeel
 * 
 */
public class Kmeans implements Clusterer {
    /**
     * The number of clusters.
     */
    private int numberOfClusters = -1;

    /**
     * The number of iterations the algorithm should make. If this value is
     * Integer.INFINITY, then the algorithm runs until the centroids no longer
     * change.
     * 
     */
    private int numberOfIterations = -1;

    /**
     * Random generator for this clusterer.
     */
    private Random rg;

    /**
     * The distance measure used in the algorithm, defaults to Euclidean
     * distance.
     */
    private DistanceMeasure dm;

    /**
     * The centroids of the different clusters.
     */
    private Instance[] centroids;
    
    /*
     * This is added by me on 4_1_16
     * Boolean to indicate whether to use a uniformly distributed initial centroids (true) or
     * use the default random initial centroid option (false)
     */
    private boolean uniform = false;
    
    /**
     * Constuct a default K-means clusterer with 100 iterations, 4 clusters, a
     * default random generator and using the Euclidean distance.
     */
    public Kmeans() {
        this(4);
    }

    /**
     * Constuct a default K-means clusterer with the specified number of
     * clusters, 100 iterations, a default random generator and using the
     * Euclidean distance.
     * 
     * @param k the number of clusters to create
     */
    public Kmeans(int k) {
        this(k, 100);
    }

    /**
     * Create a new Simple K-means clusterer with the given number of clusters
     * and iterations. The internal random generator is a new one based upon the
     * current system time. For the distance we use the Euclidean n-space
     * distance.
     * 
     * @param clusters
     *            the number of clusters
     * @param iterations
     *            the number of iterations
     */
    public Kmeans(int clusters, int iterations) {
        this(clusters, iterations, new EuclideanDistance());
    }

    /**
     * Create a new K-means clusterer with the given number of clusters and
     * iterations. Also the Random Generator for the clusterer is given as
     * parameter.
     * 
     * @param clusters
     *            the number of clustesr
     * @param iterations
     *            the number of iterations
     * 
     * @param dm
     *            the distance measure to use
     */
    public Kmeans(int clusters, int iterations, DistanceMeasure dm) {
        this.numberOfClusters = clusters;
        this.numberOfIterations = iterations;
        this.dm = dm;
        rg = new Random(System.currentTimeMillis());
    }
    
    /**
     * New method added 4-1-16. To set the uniform variable to true
     */
    public void setUniformInitialCentroids(){
    	uniform = true;
    }

    /**
     * Execute the KMeans clustering algorithm on the data set that is provided.
     * 
     * @param data data set to cluster
     * @param clusters as an array of Datasets. Each Dataset represents a cluster.
     */
    public Dataset[] cluster(Dataset data) {
        if (data.size() == 0)
            throw new RuntimeException("The dataset should not be empty");
        if (numberOfClusters == 0)
            throw new RuntimeException("There should be at least one cluster");
        // Place K points into the space represented by the objects that are
        // being clustered. These points represent the initial group of
        // centroids.
        // DatasetTools.
        Instance min = DatasetTools.minAttributes(data);
        Instance max = DatasetTools.maxAttributes(data);
        this.centroids = new Instance[numberOfClusters];
        int instanceLength = data.instance(0).noAttributes();
        if (!uniform){
        	for (int j = 0; j < numberOfClusters; j++) {
//            	double[] randomInstance = new double[instanceLength];
//            	for (int i = 0; i < instanceLength; i++) {
//              	  double dist = Math.abs(max.value(i) - min.value(i));
//              	  randomInstance[i] = (float) (min.value(i) + rg.nextDouble() * dist);
//
//          	  }
        		double[] randomInstance = DatasetTools.getRandomInstance(data, rg);
            	this.centroids[j] = new DenseInstance(randomInstance);
        	}
        }
        	
        else if(uniform){
        	DenseInstance average = (DenseInstance) DatasetTools.average(data);
        	DenseInstance stddev = (DenseInstance) DatasetTools.standardDeviation(data, average);
        	ArrayList<Double> means = (ArrayList<Double>) average.values();
        	ArrayList<Double> stddevs = (ArrayList<Double>) stddev.values();
        	
        	
        		double[][] newMeans = new double[numberOfClusters][means.size()];
        		for (int i = 0;i < means.size();i++){
        			double meanLower = means.get(i) - (3.0*stddevs.get(i));
    				double meanUpper = means.get(i) + (3.0*stddevs.get(i));
    				double meanStep = (meanUpper - meanLower) / ((double)numberOfClusters - 1);
    				for (int a = 0;a < numberOfClusters;a++){
    					newMeans[a][i] = meanLower + (a * meanStep);
    				}
        		}
        		for (int i = 0; i < newMeans.length;i++){
        			this.centroids[i] = new DenseInstance(newMeans[i]);
        		}
        	
        }
        int iterationCount = 0;
        boolean centroidsChanged = true;
        boolean randomCentroids = true;
        while (randomCentroids || (iterationCount < this.numberOfIterations && centroidsChanged)) {
        	//System.out.println("Kmeans iteration:\t"+iterationCount);
        	//System.out.println("RandomCentroids:\t"+randomCentroids);
        	//System.out.println("CentroidsChanged:\t"+centroidsChanged);
            iterationCount++;
            // Assign each object to the group that has the closest centroid.
            int[] assignment = new int[data.size()];
            for (int i = 0; i < data.size(); i++) {
                int tmpCluster = 0;
                double minDistance = dm.measure(centroids[0], data.instance(i));
                for (int j = 1; j < centroids.length; j++) {
                    double dist = dm.measure(centroids[j], data.instance(i));
                    if (dm.compare(dist, minDistance)) {
                        minDistance = dist;
                        tmpCluster = j;
                    }
                }
                assignment[i] = tmpCluster;

            }
            // When all objects have been assigned, recalculate the positions of
            // the K centroids and start over.
            // The new position of the centroid is the weighted center of the
            // current cluster.
            double[][] sumPosition = new double[this.numberOfClusters][instanceLength];
            int[] countPosition = new int[this.numberOfClusters];
            for (int i = 0; i < data.size(); i++) {
                Instance in = data.instance(i);
                for (int j = 0; j < instanceLength; j++) {

                    sumPosition[assignment[i]][j] += in.value(j);

                }
                countPosition[assignment[i]]++;
            }
            centroidsChanged = false;
            randomCentroids = false;
            for (int i = 0; i < this.numberOfClusters; i++) {
                if (countPosition[i] > 0) {
                    double[] tmp = new double[instanceLength];
                    for (int j = 0; j < instanceLength; j++) {
                        tmp[j] = (float) sumPosition[i][j] / countPosition[i];
                    }
                    Instance newCentroid = new DenseInstance(tmp);
                    if (dm.measure(newCentroid, centroids[i]) > 0.0001) {
                        centroidsChanged = true;
                        centroids[i] = newCentroid;
                    }
                } else {
                    double[] randomInstance = new double[instanceLength];
                    for (int j = 0; j < instanceLength; j++) {
                        double dist = Math.abs(max.value(j) - min.value(j));
                        randomInstance[j] = (float) (min.value(j) + rg.nextDouble() * dist);

                    }
                    randomCentroids = true;
                    this.centroids[i] = new DenseInstance(randomInstance);
                }

            }

        }
        Dataset[] output = new Dataset[centroids.length];
        for (int i = 0; i < centroids.length; i++)
            output[i] = new DefaultDataset();
        for (int i = 0; i < data.size(); i++) {
            int tmpCluster = 0;
            double minDistance = dm.measure(centroids[0], data.instance(i));
            for (int j = 0; j < centroids.length; j++) {
                double dist = dm.measure(centroids[j], data.instance(i));
                if (dm.compare(dist, minDistance)) {
                    minDistance = dist;
                    tmpCluster = j;
                }
            }
            output[tmpCluster].add(data.instance(i));

        }
        return output;
    }

}