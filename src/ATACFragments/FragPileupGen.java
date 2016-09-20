package ATACFragments;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import Node.TagNode;



public class FragPileupGen {
	private ArrayList<TagNode> frags;
	private ArrayList<TagNode> genome;
	private ArrayList<double[]> tracks;
	private int minMapQ;
	//private boolean keepDups;
	AbstractRealDistribution shortDist;
	AbstractRealDistribution monoDist;
	AbstractRealDistribution diDist;
	AbstractRealDistribution triDist;
	AbstractRealDistribution quadDist;
	
	
	public FragPileupGen(ArrayList<TagNode> f,ArrayList<TagNode> g,double[] mode, double[] means,double[] lamda) throws FileNotFoundException{
		frags = f;
		genome = g;
		shortDist = getDist(mode[0],means[0],lamda[0]);
		monoDist = getDist(mode[1],means[1],lamda[1]);
		diDist = getDist(mode[2],means[2],lamda[2]);
		triDist = getDist(mode[3],means[3],lamda[3]);
		HashMap<String,ArrayList<TagNode>> map = buildMap();
		printTracks2(map);
	}
	public FragPileupGen(File input,File index,ArrayList<TagNode> g,double[] mode, double[] means,double[] lamda,
			int q) throws FileNotFoundException{
		genome = g;
		shortDist = getDist(mode[0],means[0],lamda[0]);
		monoDist = getDist(mode[1],means[1],lamda[1]);
		diDist = getDist(mode[2],means[2],lamda[2]);
		triDist = getDist(mode[3],means[3],lamda[3]);
		minMapQ=q;
		
		printTracks3(input,index);
	}
	public ArrayList<double[]> getTracks(){return tracks;}
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
	public ArrayList<double[]> getAverageTracks(){
		ArrayList<double[]> newTracks = new ArrayList<double[]>();
		//System.out.println(tracks.size());
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
	private void printTracks3(File input,File index) throws FileNotFoundException{
		tracks = new ArrayList<double[]>();
		SAMFileReader reader = new SAMFileReader(input,index);
		
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int bedStart = genome.get(i).getStart();
			int bedStop = genome.get(i).getStop();
			//String file = chr+"_"+bedStart+"_"+bedStop+".csv";
			//PrintStream out = new PrintStream(file);
			HashMap<Integer,double[]> pileup = new HashMap<Integer,double[]>();
			for (int a = bedStart;a < bedStop;a++){
				if (!pileup.containsKey(a)){
					double[] t = new double[4];
					pileup.put(a, t);
				}
			}
			int newBedStart = bedStart;
			if (bedStart - 1000 > 0){
				newBedStart = bedStart - 1000;
			}
			int newBedStop = bedStop + 1000;
			//System.out.println(chr+newBedStart+newBedStop);
			CloseableIterator<SAMRecord> iter = reader.query(chr, bedStart, bedStop, false);
			while (iter.hasNext()){
				SAMRecord record = iter.next();
				//System.out.println("within sam iter");
				
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
								//System.out.println(t[0]+"\t"+t[1]+"\t"+t[2]+"\t"+t[3]);
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
								//System.out.println(t[0]+"\t"+t[1]+"\t"+t[2]+"\t"+t[3]);
							}
						}
					}
					
						
				}
			}
			iter.close();
			for (int x = bedStart;x < bedStop;x++){
				tracks.add(pileup.get(x));
				/*
				
				*/
			}
			/*
			for (int a = bedStart;a < bedStop;a++){
				double[] t = pileup.get(a);
				for (int x = 0 ; x < t.length;x++){
					if (x < t.length-1){
						System.out.print((t[x])+",");
					}
					else{
						System.out.println(t[x]);
					}
				}
			}
			*/
		}
		reader.close();
	}
	private AbstractRealDistribution getDist(double p,double m, double l){
		
		if (p == 2){
			return new NormalDistribution(m,l);
		}
		if (p == 0.5 || p == 3){
			return new ExponentialDistribution(m);
		}
		
		else{
			return null;
		}
		
	}
	private double greatest(double one,double two,double three,double four){
		if (one > two && one > three && one > four){
			return one;
		}
		if (two > one && two > three && two > four){
			return two;
		}
		if (three > one && three > two && three > four){
			return three;
		}
		if (four > one && four > two && four > three)
			return four;
		else
			return -1.0;
					
	}
	public void printTracks2(HashMap<String,ArrayList<TagNode>> map) throws FileNotFoundException{
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int bedStart = genome.get(i).getStart();
			int bedStop = genome.get(i).getStop();
			ArrayList<TagNode> temp = map.get(chr);
			if (temp == null){i++;}
			String file = chr+"_"+bedStart+"_"+bedStop+".csv";
			PrintStream out = new PrintStream(file);
			HashMap<Integer,double[]> pileup = new HashMap<Integer,double[]>();
			for (int a = bedStart;a < bedStop;a++){
				if (!pileup.containsKey(a)){
					double[] t = new double[4];
					pileup.put(a, t);
				}
			}
			for (int a = 0; a < temp.size();a++){
				int start = temp.get(a).getStart();
				int stop = temp.get(a).getStop();
				int length = temp.get(a).getLength();
				
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
			for (int a = bedStart;a < bedStop;a++){
				double[] t = pileup.get(a);
				for (int x = 0 ; x < t.length;x++){
					if (x < t.length-1){
						out.print((t[x])+",");
					}
					else{
						out.println(t[x]);
					}
				}
			}
		}
		
		
	}
	public void printTracks(HashMap<String,ArrayList<TagNode>> map) throws FileNotFoundException{
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int bedStart = genome.get(i).getStart();
			int bedStop = genome.get(i).getStop();
			ArrayList<TagNode> temp = map.get(chr);
			if (temp == null){i++;}
			String file = chr+"_"+bedStart+"_"+bedStop+".csv";
			PrintStream out = new PrintStream(file);
			double[][] pileup = new double[bedStop - bedStart][4];
			for (int a = 0; a < temp.size();a++){
				int start = temp.get(a).getStart();
				int stop = temp.get(a).getStop();
				//if ((bedStart <= start && bedStop >= start) || (bedStop >= stop && bedStart <= stop)){
					int begin = start;
					int end = stop;
					/*if (start < bedStart){
						begin = bedStart;
					}
					if (stop > bedStop){
						end = bedStop;
					}*/
					
					int length = temp.get(a).getLength();
					double denom = shortDist.density(length)+monoDist.density(length)+
							diDist.density(length)+triDist.density(length);
					double sh = shortDist.density(length); // denom;
					double mono = monoDist.density(length); // denom;
					double di = diDist.density(length); // denom;
					double tri = triDist.density(length); // denom;
					//double greatest = greatest(sh,mono,di,tri);
					for (int x = begin;x < end;x++ ){
						pileup[x][0] += sh;
						pileup[x][1] += mono;
						pileup[x][2] += di;
						pileup[x][3] += tri;
					
					}
			}
			for (int a = 0 ; a < pileup.length;a++){
				double d = 0.0;
				for (int x = 0;x < pileup[a].length;x++){
					 d += pileup[a][x];
				}
				for (int x = 0 ; x < pileup[a].length;x++){
					if (x < pileup[a].length-1){
						out.print((pileup[a][x])+",");
					}
					else{
						out.println(pileup[a][x]);
					}
				}
			}
			pileup = null;
			}
			
		}
	
	
	private HashMap<String,ArrayList<TagNode>> buildMap(){
		HashMap<String,ArrayList<TagNode>> map = new HashMap<String,ArrayList<TagNode>>();
		
		for (int i = 0; i < frags.size();i++){
			String chr = frags.get(i).getChrom();
			if (!map.containsKey(chr)){
				ArrayList<TagNode> temp = new ArrayList<TagNode>();
				temp.add(frags.get(i));
				map.put(chr, temp);
			}
			else{
				ArrayList<TagNode> temp = map.get(chr);
				temp.add(frags.get(i));
				map.put(chr, temp);
			}
			//frags.remove(i);
			//i--;
		}
		frags = null;
		return map;
	}
	
	private void run(){
		
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			
		}
	
	}
	
	//for testing
	public static void main(String[] args){
		NormalDistribution triDist = new NormalDistribution(600,33);
		System.out.println(triDist.density(1871));
	}
	
}

