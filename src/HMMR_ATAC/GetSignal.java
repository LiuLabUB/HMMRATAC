package HMMR_ATAC;

import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math.stat.descriptive.rank.Max;
import org.apache.commons.math.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

public class GetSignal {

	
	private  String wig = null;
	private  int start;
	private  int stop;
	private String chr;
	
	private double score;
	private double median;
	private double max;
	
	public GetSignal(String w,String c, int b,int e) throws IOException{
		wig = w;
		chr = c;
		start = b;
		stop = e;
		find();
	}
	public double getScore(){return score;}
	public double getMedian(){return median;}
	public double getMax(){return max;}
	
	private void find() throws IOException{
		BBFileReader wigReader = new BBFileReader(wig);
		BigWigIterator iter = wigReader.getBigWigIterator(chr,start,chr,stop,false);
		Max m = new Max();
		Mean mu = new Mean();
		ArrayList<Double> values = new ArrayList<Double>();
		while (iter.hasNext()){
			WigItem item = iter.next();
			double value = item.getWigValue();
			m.increment(value);
			for (int i = item.getStartBase();i < item.getEndBase();i++){
				values.add(value);
				mu.increment(value);
			}
			
		}
		if (values.size()>0){
			score = mu.getResult();
			max = m.getResult();
			median = values.get(values.size()/2);
		}
		else{
			score=0;max=0;median=0;
		}
	}
	
	
}
