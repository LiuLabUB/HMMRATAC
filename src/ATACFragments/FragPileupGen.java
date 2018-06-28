package ATACFragments;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import Node.TagNode;



public class FragPileupGen {
	private ArrayList<TagNode> genome;
	private ArrayList<double[]> tracks;
	private int minMapQ;
	AbstractRealDistribution shortDist;
	AbstractRealDistribution monoDist;
	AbstractRealDistribution diDist;
	AbstractRealDistribution triDist;
	AbstractRealDistribution quadDist;
	
	/**
	 * Constructor for creating a new FragPileupGen object, which builds the data tracks
	 * @param input a BAM file containing ATAC-seq reads
	 * @param index a BAM index file for the BAM file
	 * @param g an ArrayList of TagNodes representing the genomic regions over which the data will be created
	 * @param mode an array of doubles representing the distributions to use for building the data tracks
	 * @param means an array of doubles representing the means of the distributions
	 * @param lamda an array of doubles representing the standard deviations of the distributions
	 * @param q an integer representing the minimum mapping quality score of the reads to use for data track generation
	 * @throws FileNotFoundException
	 */
	public FragPileupGen(File input,File index,ArrayList<TagNode> g,double[] mode, double[] means,double[] lamda,
			int q) throws FileNotFoundException{
		genome = g;
		shortDist = getDist(mode[0],means[0],lamda[0]);
		monoDist = getDist(mode[1],means[1],lamda[1]);
		diDist = getDist(mode[2],means[2],lamda[2]);
		triDist = getDist(mode[3],means[3],lamda[3]);
		minMapQ=q;
		
		buildTracks(input,index);
	}
	
	/**
	 * Access the completed tracks.
	 * @return an ArrayList containing the track data as an array of doubles.
	 */
	public ArrayList<double[]> getTracks(){return tracks;}
	/**
	 * Reverse the tracks. Only used as a test.
	 * @param ori an ArrayList containing the data to be transformed as an array of doubles.
	 * @return a new ArrayList of arrays of doubles containing the transformed data.
	 */ 
	public ArrayList<double[]> reverseTracks(ArrayList<double[]> ori){
		ArrayList<double[]> reversed = new ArrayList<double[]>();
		for(int i = ori.size()-1;i >=0;i--){
			reversed.add(ori.get(i));
		}
		return reversed;
	}
	
	/**
	 * Perform a square root transformation of the tracks.
	 * @param ori an ArrayList containing the data to be transformed as an array of doubles.
	 * @return a new ArrayList of arrays of doubles containing the transformed data.
	 */
	public ArrayList<double[]> transformTracks(ArrayList<double[]> ori){
		ArrayList<double[]> trans = new ArrayList<double[]>();
		for (int i = 0;i < ori.size();i++){
			double[] temp = new double[4];
			temp[0] = Math.sqrt(ori.get(i)[0]);
			temp[1] = Math.sqrt(ori.get(i)[1]);
			temp[2] = Math.sqrt(ori.get(i)[2]);
			temp[3] = Math.sqrt(ori.get(i)[3]);
			trans.add(temp);
		}
		return trans;
	}
	/**
	 * Averages the data in the tracks across a 10bp window
	 * @return a new ArrayList containing the averaged data
	 */
	public ArrayList<double[]> getAverageTracks(){
		ArrayList<double[]> newTracks = new ArrayList<double[]>();
		for (int i = 0;i < tracks.size();i+=10){
			double[] temp = new double[4];
			double counter = 0.0;
			int end = i + 10;
			if (end > tracks.size()){
				end = (tracks.size() - end) + end;
			}
			for (int a = i; a < end;a++){
				counter+=1.0;
				temp[0]+=tracks.get(a)[0];
				temp[1]+=tracks.get(a)[1];
				temp[2]+=tracks.get(a)[2];
				temp[3]+=tracks.get(a)[3];
			}
			for (int a = 0;a < temp.length;a++){
				temp[a]/=counter;
			}
			newTracks.add(temp);
		
			
			
			
		}
		
		return newTracks;
		
	}
	/**
	 * Builds the data tracks.
	 * @param input a BAM file containing the ATAC-seq paired end reads.
	 * @param index a BAM index file containing the index for the BAM file.
	 * @throws FileNotFoundException
	 */
	private void buildTracks(File input,File index) throws FileNotFoundException{
		tracks = new ArrayList<double[]>();
		SAMFileReader reader = new SAMFileReader(input,index);
		
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int bedStart = genome.get(i).getStart();
			int bedStop = genome.get(i).getStop();
			
			HashMap<Integer,double[]> pileup = new HashMap<Integer,double[]>();
			for (int a = bedStart;a < bedStop;a++){
				if (!pileup.containsKey(a)){
					double[] t = new double[4];
					pileup.put(a, t);
				}
			}
			
			CloseableIterator<SAMRecord> iter = reader.query(chr, bedStart, bedStop, false);
			while (iter.hasNext()){
				SAMRecord record = null;
				try{
					record = iter.next();
				}
				catch(SAMFormatException ex){
					System.out.println("SAM Record is problematic. Has mapQ != 0 for unmapped read. Will continue anyway");
				}
				if(record != null){
				if(!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()
						&& record.getMappingQuality()>=minMapQ && !record.getDuplicateReadFlag()) {
					int start;
					int stop;
					if(record.getInferredInsertSize() > 0) {
						
						start = record.getAlignmentStart();
						
						stop = record.getAlignmentStart() + record.getInferredInsertSize() - 1;
						int length = stop - start;
						double sh = shortDist.density(length); 
						double mono = monoDist.density(length); 
						double di = diDist.density(length);
						double tri = triDist.density(length); 
						
						for (int x = start; x < stop;x++){
							if(pileup.containsKey(x)){
								double[] t = pileup.get(x);
								t[0] += sh;
								t[1] += mono;
								t[2] += di;
								t[3] += tri;
								pileup.put(x, t);
								
							}
						}
						
								
					} else if (record.getInferredInsertSize() < 0 ) {
						start = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
						stop = record.getAlignmentEnd();	
						int length = stop - start;
						double sh = shortDist.density(length); 
						double mono = monoDist.density(length); 
						double di = diDist.density(length);
						double tri = triDist.density(length); 
						
						for (int x = start; x < stop;x++){
							if(pileup.containsKey(x)){
								double[] t = pileup.get(x);
								t[0] += sh;
								t[1] += mono;
								t[2] += di;
								t[3] += tri;
								pileup.put(x, t);
							}
						}
					}
					
						
				}
			}
			}
			iter.close();
			for (int x = bedStart;x < bedStop;x++){
				tracks.add(pileup.get(x));
				//System.out.println("tracks\t"+chr+"\t"+x+"\t"+(x+1)+"\t"+pileup.get(x));//Added on 9/29/17 to print out actual tracks for paper. REMOVE THIS LINE WHEN DONE!!!!
				
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
	private AbstractRealDistribution getDist(double p,double m, double l){
		if (p == 1){
			//return new LaplaceDistribution(m,l);
			return new ExponentialDistribution(m);
		}
		if (p == 2){
			return new NormalDistribution(m,l);
		}
		if (p == 0.5 || p == 3){
			return new ExponentialDistribution(m);
		}
		if (p == 0){
			//return new ModifiedLaplaceDistribution(m,l);
			return new ExponentialDistribution(m);
		}
		else{
			return null;
		}
		
	}
	
	
	
}
