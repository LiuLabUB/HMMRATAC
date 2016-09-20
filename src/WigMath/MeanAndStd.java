package WigMath;

import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class MeanAndStd {
	
	private String wigFile;
	private double mean;
	private double std;

	public MeanAndStd(String w) throws IOException{
		wigFile = w;
		SetMeanAndStd();
	}
	public MeanAndStd(String w, TagNode node) throws IOException{
		wigFile = w;
		SetMeanAndStd2(node);
	}
	private void SetMeanAndStd2(TagNode node) throws IOException{
		BBFileReader wigReader = new BBFileReader(wigFile);
		String chrom = node.getChrom();
		int begin = node.getStart();
		int end = node.getStop();
		BigWigIterator iter = wigReader.getBigWigIterator(chrom, begin, chrom, end, false);
		Mean mu = new Mean();
		StandardDeviation dev = new StandardDeviation();
		while(iter.hasNext()){
			WigItem item = iter.next();
			String chr = item.getChromosome();
			int start = item.getStartBase();
			int stop = item.getEndBase();
			double value = item.getWigValue();
			for (int i = start; i < stop;i++){
				mu.increment(value);
				dev.increment(value);
			}
			
		}
		mean = mu.getResult();
		std = dev.getResult();
		
	}
	public double getMean(){return mean;}
	public double getStd(){return std;}
	private void SetMeanAndStd() throws IOException{
		BBFileReader wigReader = new BBFileReader(wigFile);
		BigWigIterator iter = wigReader.getBigWigIterator();
		Mean mu = new Mean();
		StandardDeviation dev = new StandardDeviation();
		while(iter.hasNext()){
			WigItem item = iter.next();
			String chr = item.getChromosome();
			int start = item.getStartBase();
			int stop = item.getEndBase();
			double value = item.getWigValue();
			for (int i = start; i < stop;i++){
				mu.increment(value);
				dev.increment(value);
			}
			
		}
		mean = mu.getResult();
		std = dev.getResult();
	}
}
