package WigMath;

import java.io.IOException;
import java.util.ArrayList;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class FoldChange {

	
	
	private String w = null;
	private double u;
	private double l;
	private ArrayList<TagNode> g;
	private ArrayList<TagNode> result;
	private double mean;
	private double std;
	/**
	 * Constructor for creating new FoldChange object and setting the data
	 * @param W a String representing the name of a bigwig file
	 * @param U a double representing the upper limit of the fold change range
	 * @param L a double representing the lower limit of the fold change range
	 * @param G an ArrayList of TagNode representing the genome
	 * @throws IOException
	 */
	public FoldChange(String W,double U,double L,ArrayList<TagNode> G) throws IOException{
		w = W;
		u = U;
		l = L;
		g = G;
		set();
	}
	/**
	 * Access the genomic mean
	 * @return a double representing the genomic mean
	 */
	public double getMean(){return mean;}
	/**
	 * Access the genomic standard deviation
	 * @return a double representing the genomic standard deviation
	 */
	public double getSTD(){return std;}
	/**
	 * Access the regions where the fold change is between the upper and lower limits
	 * @return an ArrayList of TagNode representing regions that meet the criteria 
	 */
	public ArrayList<TagNode> getResults(){return result;}
	/**
	 * Set the data
	 * @throws IOException
	 */
	private void set() throws IOException{
		MeanAndStd results = new MeanAndStd(w);
		mean = results.getMean();
		std = results.getStd();
		result = new ArrayList<TagNode>();
		BBFileReader wigReader = new BBFileReader(w);
		for (int i = 0;i < g.size();i++){
			String chrom = g.get(i).getChrom();
			int begin = g.get(i).getStart();
			int end = g.get(i).getStop();
			BigWigIterator iter = wigReader.getBigWigIterator(chrom,begin,chrom,end,false);
			while(iter.hasNext()){
				WigItem item = iter.next();
				int start = item.getStartBase();
				int stop = item.getEndBase();
				double value = item.getWigValue();
				if ((value/mean) >= l && (value/mean) <= u){
					TagNode temp = new TagNode(chrom,start,stop,(value/mean));
					result.add(temp);
				}
				
			}
		}
		wigReader=null;
	}
	
}
