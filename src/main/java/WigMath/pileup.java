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

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.LaplaceDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import FormatConverters.PileupToBedGraph;
import Node.PileupNode2;
import Node.TagNode;

public class pileup {

	private static ArrayList<TagNode> intervals;
	private static double start;
	private static File input;
	private static File index;
	private static int minMapQ;
	private static boolean rmDup;
	
	AbstractRealDistribution shortDist;
	AbstractRealDistribution monoDist;
	AbstractRealDistribution diDist;
	AbstractRealDistribution triDist;
	private static double shortStart;
	private static double monoStart;
	private static double diStart;
	private static double triStart;
	
	private static ArrayList<TagNode> bdg;
	private static ArrayList<TagNode> shortBG;
	private static ArrayList<TagNode> monoBG;
	private static ArrayList<TagNode> diBG;
	private static ArrayList<TagNode> triBG;
	private static double cpmScale=0;
	
	public pileup(ArrayList<TagNode> t, int s, File b, File i,int q,boolean r){
		
		intervals = t;
		start = s;
		input = b;
		index = i;
		minMapQ = q;
		rmDup=r;
		build();
	}
	public pileup(ArrayList<TagNode> t, int s, File b, File i,int q,
			double[] mode, double[] means,double[] lamda){
		intervals = t;
		start = s;
		input = b;
		index = i;
		minMapQ = q;
		shortStart=s;
		monoStart=s;
		diStart=s;
		triStart=s;
		shortDist = getDist(mode[0],means[0],lamda[0]);
		monoDist = getDist(mode[1],means[1],lamda[1]);
		diDist = getDist(mode[2],means[2],lamda[2]);
		triDist = getDist(mode[3],means[3],lamda[3]);
		buildHMMRTrack();
		
	}
	public double getCPMScale(){return cpmScale;}
	public ArrayList<TagNode> getBedGraph(){return bdg;}
	public ArrayList<TagNode> getShortBedGraph(){return shortBG;}
	public ArrayList<TagNode> getMonoBedGraph(){return monoBG;}
	public ArrayList<TagNode> getDiBedGraph(){return diBG;}
	public ArrayList<TagNode> getTriBedGraph(){return triBG;}
	
	public void buildHMMRTrack(){
		bdg = new ArrayList<TagNode>();
		shortBG = new ArrayList<TagNode>();
		monoBG= new ArrayList<TagNode>();
		diBG= new ArrayList<TagNode>();
		triBG= new ArrayList<TagNode>();
		for (int i = 0; i < intervals.size();i++){
			ArrayList<double[]> temp = makeBlockComplete(intervals.get(i));
			bdg.addAll(toBedGraph(intervals.get(i),temp.get(0)));
			shortBG.addAll(toBedGraph(intervals.get(i),temp.get(1)));
			monoBG.addAll(toBedGraph(intervals.get(i),temp.get(2)));
			diBG.addAll(toBedGraph(intervals.get(i),temp.get(3)));
			triBG.addAll(toBedGraph(intervals.get(i),temp.get(4)));
		}
	}
	public ArrayList<double[]> makeBlockComplete(TagNode t){
		double[] temp = new double[t.getLength()];
		double[] tempShort = new double[t.getLength()];
		double[] tempMono = new double[t.getLength()];
		double[] tempDi = new double[t.getLength()];
		double[] tempTri = new double[t.getLength()];
		
		SAMFileReader reader = new SAMFileReader(input,index);
		CloseableIterator<SAMRecord> iter = reader.query(t.getChrom(), t.getStart(), 
				t.getStop(), false);
		while (iter.hasNext()){
			SAMRecord record = iter.next();
			if(!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()
					//&& record.getMappingQuality()>=minMapQ 
					&& !record.getDuplicateReadFlag()) {
				int readStart = record.getAlignmentStart();
				int readStop= record.getAlignmentStart() + record.getInferredInsertSize() - 1;
				if (record.getInferredInsertSize() < 0 ) {
					readStart = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
					readStop = record.getAlignmentEnd();	
				}
				int length = readStop - readStart;
				if (readStart < t.getStart()){
					readStart = t.getStart();
				}
				if (readStop < t.getStart()){
					readStop = t.getStart();
				}
				if (readStop >= t.getStop()){
					readStop = t.getStop()-1;
				}
				if (readStart >= t.getStop()){
					readStart = t.getStop()-1;
				}
				temp[readStart - t.getStart()]++;
				temp[readStop - t.getStart()]--;
				
				tempShort[readStart - t.getStart()]+=shortDist.density(length);
				tempShort[readStop - t.getStart()]-=shortDist.density(length);
				
				tempMono[readStart - t.getStart()]+=monoDist.density(length);
				tempMono[readStop - t.getStart()]-=monoDist.density(length);
				
				tempDi[readStart - t.getStart()]+=diDist.density(length);
				tempDi[readStop - t.getStart()]-=diDist.density(length);
				
				tempTri[readStart - t.getStart()]+=triDist.density(length);
				tempTri[readStop - t.getStart()]-=triDist.density(length);
				
					
			}
		}
		iter.close();
		reader.close();
		for (int i = 0;i < temp.length;i++){
			if (i == 0){
				temp[i] = temp[i]+start;
				tempShort[i] = tempShort[i]+shortStart;
				tempMono[i] = tempMono[i]+monoStart;
				tempDi[i] = tempDi[i]+diStart;
				tempTri[i] = tempTri[i]+triStart;
			} else{
				temp[i] = temp[i] + temp[i-1];
				tempShort[i] = tempShort[i] + tempShort[i-1];
				tempMono[i] = tempMono[i] + tempMono[i-1];
				tempDi[i] = tempDi[i] + tempDi[i-1];
				tempTri[i] = tempTri[i] + tempTri[i-1];
				
			}
		}
		
		start = temp[temp.length-1];
		shortStart = tempShort[tempShort.length-1];
		monoStart = tempMono[tempMono.length-1];
		diStart = tempDi[tempDi.length-1];
		triStart = tempDi[tempDi.length-1];
		
		ArrayList<double[]> result = new ArrayList<double[]>();
		result.add(temp);
		result.add(tempShort);
		result.add(tempMono);
		result.add(tempDi);
		result.add(tempTri);
		
		return result;
	}
	public void build(){
		bdg = new ArrayList<TagNode>();
		HashMap<String,ArrayList<TagNode>> temp = bedGraphMath.toMap(intervals);
		for (String chr : temp.keySet()){
			ArrayList<TagNode> temp1 = temp.get(chr);
			Collections.sort(temp1,TagNode.basepairComparator);
			for (int i = 0; i < temp1.size();i++){
				bdg.addAll(toBedGraph(temp1.get(i),makeBlock(temp1.get(i))));
			}
		}
	}
	public ArrayList<TagNode> toBedGraph(TagNode t,double[] temp){
		ArrayList<TagNode> bedGraph = new ArrayList<TagNode>();
		ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
		for (int i = 0; i < temp.length;i++){
			PileupNode2 pNode = new PileupNode2(t.getStart()+i,temp[i],t.getChrom());
			pile.add(pNode);
		}
		bedGraph.addAll(new PileupToBedGraph(pile,1).getBedGraph());
		
		return bedGraph;
	}
	
	public double[] makeBlock(TagNode t){
		double[] temp = new double[t.getLength()];
		
		SAMFileReader reader = new SAMFileReader(input,index);
		CloseableIterator<SAMRecord> iter = reader.query(t.getChrom(), t.getStart(), 
				t.getStop(), false);
		
		while (iter.hasNext()){
			SAMRecord record = null;
			try{
				record = iter.next();
			}
			catch(SAMFormatException ex){
				System.out.println("SAM Record is problematic. Has mapQ != 0 for unmapped read. Will continue anyway");
			}
			if(record != null){
				if (!record.getReadUnmappedFlag()
						&& !record.getMateUnmappedFlag()
						&& record.getMappingQuality() >= minMapQ
						&& !(record.getDuplicateReadFlag() && rmDup)
						&& Math.abs(record.getInferredInsertSize()) <= 1000
						&& record.getInferredInsertSize() != 0) {
					int readStart = record.getAlignmentStart();
					int readStop = record.getAlignmentEnd();
					//				int readStop= record.getAlignmentStart() + record.getInferredInsertSize() - 1;
					//				if (record.getInferredInsertSize() < 0 ) {
					//					readStart = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
					//					readStop = record.getAlignmentEnd();	
					//				}

					cpmScale++;
					if (readStart < t.getStart()) {
						readStart = t.getStart();
					}
					if (readStop < t.getStart()) {
						readStop = t.getStart();
					}
					if (readStop >= t.getStop()) {
						readStop = t.getStop() - 1;
					}
					if (readStart >= t.getStop()) {
						readStart = t.getStop() - 1;
					}
					temp[readStart - t.getStart()]++;
					temp[readStop - t.getStart()]--;

				}
			}
			
		}
		iter.close();
		reader.close();
		for (int i = 0;i < temp.length;i++){
			if (i == 0){
				temp[i] = temp[i];
			} else{
				temp[i] = temp[i] + temp[i-1];
				
			}
		}
		
		//start = temp[temp.length-1];
		
		return temp;
	}
	private AbstractRealDistribution getDist(double p,double m, double l){
		if (p == 1){
			return new LaplaceDistribution(m,l);
			//return new ExponentialDistribution(m);
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
