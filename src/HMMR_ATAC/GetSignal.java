package HMMR_ATAC;

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
import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class GetSignal {

	
	private  String wig = null;
	private  int start;
	private  int stop;
	private String chr;
	
	private double score;
	private double median;
	private double max;
	private ArrayList<Double> val;
	private TagNode node;
	/**
	 * Constructor
	 * @param w a String representing the bigwig file
	 * @param c a String representing the chromosome name
	 * @param b an integer representing the start of the region to score
	 * @param e an integer representing the stop of the region to score
	 * @throws IOException
	 */
	public GetSignal(String w,String c, int b,int e) throws IOException{
		wig = w;
		chr = c;
		start = b;
		stop = e;
		find();
	}
	/**
	 * Access the average score
	 * @return a double representing the average score
	 */
	public double getScore(){return score;}
	/**
	 * Access the median score
	 * @return a double representing the median score
	 */
	public double getMedian(){return median;}
	/**
	 * Access the max score
	 * @return a double representing the maximum score
	 */
	public double getMax(){return max;}
	/**
	 * Access the bigwig scores
	 * @return an ArrayList of doubles representing all the bigwig scores 
	 */
	public ArrayList<Double> getValues(){return val;}
	/**
	 * Access the position with the maximum score
	 * @return a TagNode representing the position having the maximum score
	 */
	public TagNode getMaxPos(){return node;}
	/**
	 * Find the max and the max position
	 * @throws IOException
	 */
	private void find() throws IOException{
		BBFileReader wigReader = new BBFileReader(wig);
		BigWigIterator iter = wigReader.getBigWigIterator(chr,start,chr,stop,false);
		Max m = new Max();
		Mean mu = new Mean();
		ArrayList<Double> values = new ArrayList<Double>();
		double tempMax=0;
		node  = null;
		while (iter.hasNext()){
			WigItem item = iter.next();
			double value = item.getWigValue();
			m.increment(value);
			for (int i = item.getStartBase();i < item.getEndBase();i++){
				values.add(value);
				mu.increment(value);
				if (value > tempMax){
					tempMax = value;
					node = new TagNode(chr,i,i+1);
				}
			}
			
		}
		
		val = values;
		if (values.size()>0){
			score = mu.getResult();
			max = m.getResult();
			median = values.get(values.size()/2);
		}
		else{
			score=0;max=0;median=0;
		}
		wigReader = null;
	}
	/**
	 * Find max positions and extend toward second highest position
	 * @throws IOException
	 */
	public void findAll() throws IOException{
		if (max > 0){
		ArrayList<TagNode> temp = new ArrayList<TagNode>();
		BBFileReader wigReader = new BBFileReader(wig);
		BigWigIterator iter = wigReader.getBigWigIterator(chr,start,chr,stop,false);
		
		while (iter.hasNext()){
			WigItem item = iter.next();
			double value = item.getWigValue();
			
			for (int i = item.getStartBase();i < item.getEndBase();i++){
				if (value == max){
					TagNode tempNode = new TagNode(chr,i,i+1);
					temp.add(tempNode);
				}
			}
			
		}
		if (temp.size() > 1){
			TagNode left = temp.get(0);
			TagNode right = temp.get(temp.size()-1);
			int start = ((right.getStop() - left.getStart())/2) + left.getStart();
			node = new TagNode(left.getChrom(),start,start+1);
		}
		wigReader=null;
		}
	}
	/**
	 * Find max and max position of gaussian smoothed data
	 * @param stdev an integer representing the standard deviation for smoothing
	 * @throws IOException
	 */
	public void findSmooth(int stdev) throws IOException{
		
          if (max > 0){  
            
            // Use a window size equal to +/- 3 SD's
            
            double[] filter = new double[6*stdev+1];
            double sum = 0;
            for (int i = 0; i < filter.length; i++) {
                    double x = i - 3*stdev;
                    double value = (double) Math.exp(-(x*x) / (2*stdev*stdev));
                    filter[i] = value;
                    sum += value;
            }
            // Normalize so that the filter is area-preserving (has total area = 1)
            for (int i = 0; i < filter.length; i++) {
                    filter[i] /= sum;
            }

            BBFileReader wigReader = new BBFileReader(wig);
            int paddedStart = start-(3*stdev);
            int paddedStop = stop+(3*stdev);
    		BigWigIterator iter = wigReader.getBigWigIterator(chr,paddedStart,chr,paddedStop,false);
    		double[] data = new double[paddedStop-paddedStart];
    		while (iter.hasNext()){
    			WigItem item = iter.next();
    			double value = item.getWigValue();
    			
    			for (int i = item.getStartBase();i < item.getEndBase();i++){
    				if (i >= paddedStart && i < paddedStop){
    					data[i - paddedStart] = value;
    				}
    			}
    			
    		}
    		wigReader=null;
    		
            // Convolve the data with the filter
            double[] smoothed = new double[stop-start];
            for (int i = 0; i < smoothed.length; i++) {
                    for (int j = 0; j < filter.length; j++) {
                            smoothed[i] += data[i+j] * filter[j];
                    }
            }
            
            double tempMax = 0.0;
            for (int i = 0;i < smoothed.length;i++){
            	if (smoothed[i] > tempMax){
            		tempMax = smoothed[i];
            		node = new TagNode(chr,start+i,start+i+1);
            	}
            }
            
	}
            
    }
		
	
	
}
