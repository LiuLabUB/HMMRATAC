package HMMR_ATAC;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import Node.TagNode;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

public class pullLargeLengths {
	
	private File bam;
	private File index;
	private int minQ;
	private double[] lengths;
	private ArrayList<TagNode> genome;
	private double[] weights = new double[3];
	private double[] means;
	public pullLargeLengths(File b, File i,int m,ArrayList<TagNode> g,double[] mu){
		bam = b;
		index = i;
		minQ = m;
		genome = g;
		means = mu;
		read();
		setWeights();
	}
	public double[] getLengths(){return lengths;}
	public double[] getWeights(){return weights;}
	public double[] getSampledLengths(int div){
		int len = lengths.length/div;
		double[] newLengths = new double[len];
		shuffle(lengths);
		for(int i = 0; i < len;i++){
			newLengths[i] = lengths[i];
		}
		return newLengths;
	}
	private void shuffle(double[] list){
		Random rnd = new Random();
		for(int i = list.length -1;i > 0;i--){
			int index = rnd.nextInt(i+1);
			double a = list[index];
			list[index] = list[i];
			list[i] = a;
		}
	}
	private void read(){
		int counter = 0;
		SAMFileReader reader = new SAMFileReader(bam,index);
		ArrayList<Double> temp = new ArrayList<Double>();
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int start = genome.get(i).getStart();
			int stop = genome.get(i).getStop();
			CloseableIterator<SAMRecord> iter = reader.query(chr,start,stop,false);
			while (iter.hasNext()){
				SAMRecord record = iter.next();
				if(!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()) {
					if (record.getMappingQuality() >= minQ){
						
						if (Math.abs(record.getInferredInsertSize()) > 100 && Math.abs(record.getInferredInsertSize())
								< 1000){
							counter+=1;
							temp.add((double)Math.abs(record.getInferredInsertSize()));
						}
					}
				}
			}
			iter.close();
		}
		reader.close();
		lengths = new double[counter];
		for (int i = 0;i < temp.size();i++){
			if (temp.get(i) > 100){
				lengths[i] = temp.get(i);
			}
		}
		
	}
	private void setWeights(){
		double cutone = (means[1] - means[0])/2 + means[0];
		double cuttwo = (means[2] - means[1])/2 + means[1];
		int counter1=0;
		int counter2=0;
		int counter3=0;
		for (int i = 0;i < lengths.length;i++){
			if (lengths[i] < cutone)
				counter1++;
			else if(lengths[i] >= cutone && lengths[i] <= cuttwo)
				counter2++;
			else
				counter3++;
		}
		weights[0] = (double)counter1/(double)lengths.length; 
		weights[1] = (double)counter2/(double)lengths.length;
		weights[2] = (double)counter3/(double)lengths.length;
		
	}

}
