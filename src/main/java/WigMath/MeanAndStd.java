package WigMath;
/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
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

	/**
	 * Constructor
	 * @param w a String representing the bigwig file
	 * @throws IOException
	 */
	public MeanAndStd(String w) throws IOException{
		wigFile = w;
		SetMeanAndStd();
	}
	/**
	 * Constructor
	 * @param w a String representing the bigwig file
	 * @param node a TagNode representing a region to calculate the mean and stddev
	 * @throws IOException
	 */
	public MeanAndStd(String w, TagNode node) throws IOException{
		wigFile = w;
		SetMeanAndStd2(node);
	}
	/**
	 * Set the data across a specific region
	 * @param node a TagNode representing a specific region for calculation
	 * @throws IOException
	 */
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
		wigReader=null;
	}
	/**
	 * Access the mean
	 * @return a double representing the mean
	 */
	public double getMean(){return mean;}
	/**
	 * Access the standard deviation
	 * @return a double representing the stddev
	 */
	public double getStd(){return std;}
	/**
	 * Set the data across the entire genome
	 * @throws IOException
	 */
	private void SetMeanAndStd() throws IOException{
		BBFileReader wigReader = new BBFileReader(wigFile);
		BigWigIterator iter = wigReader.getBigWigIterator();
		Mean mu = new Mean();
		StandardDeviation dev = new StandardDeviation();
		while(iter.hasNext()){
			WigItem item = iter.next();
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
