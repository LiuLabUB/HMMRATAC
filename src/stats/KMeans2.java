package stats;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.stat.correlation.StorelessCovariance;

import Node.ATACClusterNode;
import Node.MatrixNodeForKMeans;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.Instance;


public class KMeans2 {
	
	private Dataset data;
	private int k;
	private int numIter;
	private Dataset[] clusteredData;
	private KMeans kmeans;
	
	private ArrayList<ATACClusterNode> clusterList;
	private ArrayList<double[][]> covMat;
	
	
	public KMeans2(Dataset DATA){
		data = DATA;
		kmeans = new KMeans();
	}
	public KMeans2(Dataset DATA,int K){
		data = DATA;
		k=K;
		kmeans = new KMeans(k);
	}
	public KMeans2(Dataset DATA,int K,int iter){
		data = DATA;
		k=K;
		numIter = iter;
		kmeans = new KMeans(k,numIter);
		//kmeans.setUniformInitialCentroids();
	}
	
	public Dataset[] cluster(){
		clusteredData = kmeans.cluster(data);
		
		return clusteredData;
		
	}
	public ArrayList<ATACClusterNode> getClusterList(){
		return clusterList;
	}
	public ArrayList<double[][]> getCovMat(){
		return covMat;
	}
	
	public void makeClusteredList2Signals(Dataset[] clusters,HashMap<Integer,MatrixNodeForKMeans> map){
		clusterList = new ArrayList<ATACClusterNode>();
		ATACClusterNode temp = null;
		covMat = new ArrayList<double[][]>();
		for (int i = 0; i < clusters.length;i++){
			Dataset cluster = clusters[i];
			StorelessCovariance cov = new StorelessCovariance(2);
			
			for (int x = 0; x< cluster.size();x++){
				Instance ins = cluster.get(x);
				MatrixNodeForKMeans node = map.get(ins.getID());
				temp = new ATACClusterNode(node.getChrom(),node.getPos(),node.getEnrich1(),node.getEnrich2(),
						node.getEnrich3(),node.getIndex(),i);
				clusterList.add(temp);
				double[] row1 = new double[2];
				row1[0] = node.getEnrich1();
				row1[1] = node.getEnrich2();
				cov.increment(row1);
			}
			double[][] covM = cov.getCovarianceMatrix().getData();
			covMat.add(covM);
		}
		clusters = null; map = null;
		Collections.sort(clusterList,ATACClusterNode.positionComparator);
	}

	public void makeClusteredList3Signals(Dataset[] clusters,HashMap<Integer,MatrixNodeForKMeans> map){
		clusterList = new ArrayList<ATACClusterNode>();
		ATACClusterNode temp = null;
		covMat = new ArrayList<double[][]>();
		for (int i = 0; i < clusters.length;i++){
			Dataset cluster = clusters[i];
			StorelessCovariance cov = new StorelessCovariance(3);
			
			for (int x = 0; x< cluster.size();x++){
				Instance ins = cluster.get(x);
				MatrixNodeForKMeans node = map.get(ins.getID());
				temp = new ATACClusterNode(node.getChrom(),node.getPos(),node.getEnrich1(),node.getEnrich2(),
						node.getEnrich3(),node.getIndex(),i);
				clusterList.add(temp);
				double[] row1 = new double[3];
				row1[0] = node.getEnrich1();
				row1[1] = node.getEnrich2();
				row1[2] = node.getEnrich3();
				cov.increment(row1);
			}
			double[][] covM = cov.getCovarianceMatrix().getData();
			covMat.add(covM);
		}
		clusters = null; map = null;
		Collections.sort(clusterList,ATACClusterNode.positionComparator);
	}
}
