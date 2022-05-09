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

import Node.TagNode;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;


public class pullLargeLengths {
	
	private final File bam;
	private final File index;
	private final int minQ;
	private double[] lengths;
	private final ArrayList<TagNode> genome;
	private final double[] weights = new double[3];
	private final double[] means;
	
	/**
	 * Constructor for creating pullLargeLengths object, reading the data and setting weights
	 *
	 * @param bam  		a File representing the BAM file of reads
	 * @param index  	a File representing the BAM index file of reads
	 * @param minQ  	an integer representing the minimum mapping quality of reads to pull
	 * @param genome	an ArrayList of TagNode representing the genome
	 * @param means	 	an Array of doubles representing the means of the fragment length distributions
	 */
	public pullLargeLengths(File bam, File index, int minQ, ArrayList<TagNode> genome, double[] means) {
		this.bam = bam;
		this.index = index;
		this.minQ = minQ;
		this.genome = genome;
		this.means = means;
		read();
		setWeights();
	}
	
	/**
	 * Access the lengths
	 *
	 * @return an Array of doubles representing the lengths of the fragments read
	 */
	public double[] getLengths() {
		return lengths;
	}
	
	/**
	 * Access the weights
	 *
	 * @return an Array of doubles representing the proportion of fragments in each length category
	 */
	public double[] getWeights() {
		return weights;
	}
	
	/**
	 * Access the down sampled weights
	 *
	 * @param div an integer representing the factor to divide the number of lengths by to downsample
	 * @return an Array of doubles representing the downsampled lengths
	 */
	public double[] getSampledLengths(int div) {
		int len = lengths.length / div;
		double[] newLengths = new double[len];
		shuffle(lengths);
		System.arraycopy(lengths, 0, newLengths, 0, len);
		return newLengths;
	}
	
	/**
	 * Shuffle the length data
	 *
	 * @param arr an Array of doubles representing the pulled fragment lengths
	 */
	private void shuffle(double[] arr) {
		Random rnd = new Random();
		for (int i = arr.length - 1; i > 0; i--) {
			int index = rnd.nextInt(i + 1);
			double a = arr[index];
			arr[index] = arr[i];
			arr[i] = a;
		}
	}
	
	/**
	 * Read the data and create a list of lengths
	 */
	private void read() {
		int counter = 0;
		SAMFileReader reader = new SAMFileReader(bam, index);
		ArrayList<Double> temp = new ArrayList<>();
		for (TagNode tagNode : genome) {
			String chr = tagNode.getChrom();
			int start = tagNode.getStart();
			int stop = tagNode.getStop();
			CloseableIterator<SAMRecord> iter = reader.query(chr, start, stop, false);
			while (iter.hasNext()) {
				SAMRecord record = null;
				try {
					record = iter.next();
				} catch (SAMFormatException ex) {
					System.out.println("SAM Record is problematic. Has mapQ != 0 for unmapped read. Will continue anyway");
				}
				if (record != null && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()
						&& record.getFirstOfPairFlag() && record.getMappingQuality() >= minQ
						&& Math.abs(record.getInferredInsertSize()) > 100
						&& Math.abs(record.getInferredInsertSize()) < 1000) {
					
					counter++;
					temp.add((double) Math.abs(record.getInferredInsertSize()));
				}
			}
			iter.close();
		}
		reader.close();
		lengths = new double[counter];
		for (int i = 0; i < temp.size(); i++) {
			if (temp.get(i) > 100) {
				lengths[i] = temp.get(i);
			}
		}
	}
	
	/**
	 * Set the weights ie the proportion of fragments in each length category
	 */
	private void setWeights() {
		double cut1 = (means[1] - means[0]) / 2 + means[0];
		double cut2 = (means[2] - means[1]) / 2 + means[1];
		int counter1 = 0;
		int counter2 = 0;
		int counter3 = 0;
		for (double length : lengths) {
			if (length < cut1)
				counter1++;
			else if (length >= cut1 && length <= cut2)
				counter2++;
			else
				counter3++;
		}
		weights[0] = (double) counter1 / (double) lengths.length;
		weights[1] = (double) counter2 / (double) lengths.length;
		weights[2] = (double) counter3 / (double) lengths.length;
	}
}
