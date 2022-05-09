package ATACFragments;

/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu
 * <p>
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

import JAHMMTest.Slope;
import Node.TagNode;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

public class FragPileupGen {
	private final ArrayList<TagNode> genome;
	private ArrayList<double[]> tracks;
	
	private final int minMapQ;
	AbstractRealDistribution shortDist;
	AbstractRealDistribution monoDist;
	AbstractRealDistribution diDist;
	AbstractRealDistribution triDist;
	private final boolean rmDup;
	private final double scale;
	
	/**
	 * Constructor for creating a new FragPileupGen object, which builds the data tracks
	 * @param input a BAM file containing ATAC-seq reads
	 * @param index a BAM index file for the BAM file
	 * @param genome an ArrayList of TagNodes representing the genomic regions over which the data will be created
	 * @param mode an array of doubles representing the distributions to use for building the data tracks
	 * @param means an array of doubles representing the means of the distributions
	 * @param lambda an array of doubles representing the standard deviations of the distributions
	 * @param minMapQ an integer representing the minimum mapping quality score of the reads to use for data track generation
	 */
	public FragPileupGen(File input, File index, ArrayList<TagNode> genome, double[] mode, double[] means,
						 double[] lambda, int minMapQ, boolean rmDup, double scale) {
		this.genome = genome;
		this.shortDist = getDist(mode[0], means[0], lambda[0]);
		this.monoDist = getDist(mode[1], means[1], lambda[1]);
		this.diDist = getDist(mode[2], means[2], lambda[2]);
		this.triDist = getDist(mode[3], means[3], lambda[3]);
		this.minMapQ = minMapQ;
		this.rmDup = rmDup;
		this.scale = scale;
		
		buildTracks(input, index);
	}
	
	/**
	 * Access the completed tracks.
	 * @return an ArrayList containing the track data as an array of doubles.
	 */
	public ArrayList<double[]> getTracks() {
		return tracks;
	}
	
	/**
	 * Reverse the tracks. Only used as a test.
	 * @param ori an ArrayList containing the data to be transformed as an array of doubles.
	 * @return a new ArrayList of arrays of doubles containing the transformed data.
	 */
	public ArrayList<double[]> reverseTracks(ArrayList<double[]> ori) {
		ArrayList<double[]> reversed = new ArrayList<>();
		Collections.copy(ori, reversed);
		Collections.reverse(reversed);
		return reversed;
	}
	
	/**
	 * Perform a transformation of the data and convert the raw data into the SLOPE of the raw data
	 * @param ori an ArrayList containing the data to be transformed
	 * @param window an integer specifying the window size for calculating the slope
	 * @return a new ArrayList of arrays of doubles containing the slopes of the ori data
	 */
	public ArrayList<double[]> getSlope(ArrayList<double[]> ori, int window) {
		ArrayList<double[]> trans = new ArrayList<>();
		double[] temp1 = new double[ori.size()];
		double[] temp2 = new double[ori.size()];
		double[] temp3 = new double[ori.size()];
		double[] temp4 = new double[ori.size()];
		for (int i = 0; i < ori.size(); i++) {
			temp1[i] = ori.get(i)[0];
			temp2[i] = ori.get(i)[1];
			temp3[i] = ori.get(i)[2];
			temp4[i] = ori.get(i)[3];
		}
		temp1 = Slope.build(temp1, window);
		temp2 = Slope.build(temp2, window);
		temp3 = Slope.build(temp3, window);
		temp4 = Slope.build(temp4, window);
		for (int i = 0; i < temp1.length; i++) {
			double[] temp = new double[4];
			temp[0] = temp1[i];
			temp[1] = temp2[i];
			temp[2] = temp3[i];
			temp[3] = temp4[i];
			trans.add(temp);
		}
		return trans;
	}
	
	/**
	 * Add two ArrayList of double arrays together
	 * @param first an ArrayList of double arrays with the first set of data for merging
	 * @param second an ArrayList of double arrays with the second set of data for merging
	 * @return An ArrayList of double arrays with the merged data
	 */
	
	public ArrayList<double[]> merge(ArrayList<double[]> first, ArrayList<double[]> second) {
		ArrayList<double[]> merged = new ArrayList<>();
		if (first.size() != second.size()) {
			return null;
		} else {
			for (int i = 0; i < first.size(); i++) {
				double[] temp = new double[8];
				temp[0] = first.get(i)[0];
				temp[1] = first.get(i)[1];
				temp[2] = first.get(i)[2];
				temp[3] = first.get(i)[3];
				temp[4] = second.get(i)[0];
				temp[5] = second.get(i)[1];
				temp[6] = second.get(i)[2];
				temp[7] = second.get(i)[3];
				merged.add(temp);
			}
			return merged;
		}
		
	}
	
	/**
	 * Perform a scaling of the tracks based on inputted scaling factor.
	 * @param ori an ArrayList containing the data to be transformed as an array of doubles.
	 * @return a new ArrayList of arrays of doubles containing the transformed data.
	 */
	public ArrayList<double[]> scaleTracks(ArrayList<double[]> ori) {
		ArrayList<double[]> trans = new ArrayList<>();
		for (double[] doubles : ori) {
			double[] temp = new double[4];
			temp[0] = (doubles[0] / scale);
			temp[1] = (doubles[1] / scale);
			temp[2] = (doubles[2] / scale);
			temp[3] = (doubles[3] / scale);
			trans.add(temp);
		}
		return trans;
	}
	
	/**
	 * Perform a square root transformation of the tracks.
	 * @param ori an ArrayList containing the data to be transformed as an array of doubles.
	 * @return a new ArrayList of arrays of doubles containing the transformed data.
	 */
	public ArrayList<double[]> transformTracks(ArrayList<double[]> ori) {
		ArrayList<double[]> trans = new ArrayList<>();
		for (double[] doubles : ori) {
			double[] temp = new double[4];
			temp[0] = Math.sqrt(doubles[0]);
			temp[1] = Math.sqrt(doubles[1]);
			temp[2] = Math.sqrt(doubles[2]);
			temp[3] = Math.sqrt(doubles[3]);
			trans.add(temp);
		}
		return trans;
	}
	
	/**
	 * Averages the data in the tracks across a 10bp window
	 * @return a new ArrayList containing the averaged data
	 */
	public ArrayList<double[]> getAverageTracks() {
		ArrayList<double[]> newTracks = new ArrayList<>();
		for (int i = 0; i < tracks.size(); i += 10) {
			double[] temp = new double[4];
			double counter = 0.0;
			int end = i + 10;
			if (end > tracks.size()) {
				end = (tracks.size() - end) + end;
			}
			for (int a = i; a < end; a++) {
				counter += 1.0;
				temp[0] += tracks.get(a)[0];
				temp[1] += tracks.get(a)[1];
				temp[2] += tracks.get(a)[2];
				temp[3] += tracks.get(a)[3];
			}
			for (int a = 0; a < temp.length; a++) {
				temp[a] /= counter;
			}
			newTracks.add(temp);
		}
		return newTracks;
	}
	
	/**
	 * Builds the data tracks.
	 * @param input a BAM file containing the ATAC-seq paired end reads.
	 * @param index a BAM index file containing the index for the BAM file.
	 */
	private void buildTracks(File input, File index) {
		tracks = new ArrayList<>();
		SAMFileReader reader = new SAMFileReader(input, index);
		
		for (TagNode tagNode : genome) {
			String chr = tagNode.getChrom();
			int bedStart = tagNode.getStart();
			int bedStop = tagNode.getStop();
			
			HashMap<Integer, double[]> pileup = new HashMap<>();
			for (int a = bedStart; a < bedStop; a++) {
				if (!pileup.containsKey(a)) {
					double[] t = new double[4];
					pileup.put(a, t);
				}
			}
			
			CloseableIterator<SAMRecord> iter = reader.query(chr, bedStart, bedStop, false);
			while (iter.hasNext()) {
				SAMRecord record = null;
				try {
					record = iter.next();
				} catch (SAMFormatException ex) {
					System.out.println("SAM Record is problematic. Has mapQ != 0 for unmapped read. Will continue anyway");
				}
				if (record != null) {
					if (!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()
							&& record.getMappingQuality() >= minMapQ && !(record.getDuplicateReadFlag() && rmDup)) {
						
						int start;
						int stop;
						if (record.getInferredInsertSize() > 0) {
							start = record.getAlignmentStart();
							stop = record.getAlignmentStart() + record.getInferredInsertSize() - 1;
						} else if (record.getInferredInsertSize() < 0) {
							start = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
							stop = record.getAlignmentEnd();
						} else {
							continue;
						}
						
						int length = stop - start;
						double sh = shortDist.density(length);
						double mono = monoDist.density(length);
						double di = diDist.density(length);
						double tri = triDist.density(length);
						
						double total = sh + mono + di + tri;
						
						if (total == 0) {
							total = 1;
						}
						
						for (int x = start; x < stop; x++) {
							if (pileup.containsKey(x)) {
								double[] t = pileup.get(x);
								t[0] += sh / total;
								t[1] += mono / total;
								t[2] += di / total;
								t[3] += tri / total;
								pileup.put(x, t);
							}
						}
					}
				}
			}
			iter.close();
			for (int x = bedStart; x < bedStop; x++) {
				tracks.add(pileup.get(x));
			}
		}
		reader.close();
	}
	
	/**
	 * Returns a new Distribution 
	 * @param p a double specifying which distribution to use
	 * @param m a double representing the mean for the desired distribution
	 * @param l a double representing the standard deviation, or more generally, the lambda of the desired distribution
	 * @return A new AbstractRealDistribution
	 */
	private AbstractRealDistribution getDist(double p, double m, double l) {
		if (p == 1) {
			return new ExponentialDistribution(m);
		}
		if (p == 2) {
			return new NormalDistribution(m, l);
		}
		if (p == 0.5 || p == 3) {
			return new ExponentialDistribution(m);
		}
		if (p == 0) {
			return new ExponentialDistribution(m);
		} else {
			return null;
		}
	}
}
