package WigMath;

import java.io.IOException;
import java.util.ArrayList;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class PullWigAboveCutoff {

	private ArrayList<TagNode> g;
	
	private String wig;
	private double cut;
	private ArrayList<TagNode> result;
	/**
	 * Constructor
	 * @param w a String representing the bigwig file
	 * @param c a double representing the zscored cutoff for reporting
	 * @param G an ArrayList of TagNode representing the regions to limit the search to
	 * @throws IOException
	 */
	public PullWigAboveCutoff(String w,double c, ArrayList<TagNode> G) throws IOException{
		wig = w;
		cut = c;
		g = G;
		set();
	}
	/**
	 * Access the regions that are above the cutoff
	 * @return
	 */
	public ArrayList<TagNode> getResults(){return result;}
	/**
	 * Set the data
	 * @throws IOException
	 */
	private void set() throws IOException{
		result = new ArrayList<TagNode>();
		
		for (int i = 0;i < g.size();i++){
			MeanAndStd results = new MeanAndStd(wig,g.get(i));
			double mean = results.getMean();
			double std = results.getStd();
			results=null;
			BBFileReader reader = new BBFileReader(wig);
			String chrom = g.get(i).getChrom();
			int begin = g.get(i).getStart();
			int end = g.get(i).getStop();
			BigWigIterator iter = reader.getBigWigIterator(chrom,begin,chrom,end,false);
			while (iter.hasNext()){
				WigItem item = iter.next();
				item.getChromosome();
				int start = item.getStartBase();
				int stop = item.getEndBase();
				double value = item.getWigValue();
				value = (value-mean) / std;
				if (value >= cut){
					TagNode temp = new TagNode(chrom,start,stop,value);
					result.add(temp);
				}
			
			}
		}
		
	}
	
}
