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

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import FormatConverters.PileupToBedGraph;
import Node.PileupNode2;
import Node.TagNode;

public class PileupGenerator {
	//inputs
	private File bam;
	private File index;
	private ArrayList<TagNode> genome;
	
	//outputs
	private ArrayList<TagNode> bedgraph;
	private double mean;
	private double stddev;
	
	public PileupGenerator(File b, File i,ArrayList<TagNode> g){
		bam = b;
		index = i;
		genome = g;
		bedgraph = new ArrayList<TagNode>();
		Generate();
	}
	public ArrayList<TagNode> getBedGraph(){return bedgraph;}
	public double getMean(){return mean;}
	public double getStdDev(){return stddev;}
	
	//need to change to account for fragments and not simply reads
	private void Generate(){
		SAMFileReader reader = new SAMFileReader(bam,index);
		Mean mu = new Mean();
		StandardDeviation std = new StandardDeviation();
		for (int i = 0; i < genome.size();i++){
			String chr = genome.get(i).getChrom();
			int start = genome.get(i).getStart();
			int stop = genome.get(i).getStop();
			ArrayList<PileupNode2> temp = new ArrayList<PileupNode2>();
			for (int a = start;a < stop;a++){
				//new
				int begin;
				int end;
				if (a >= 1000){
					begin = a -1000;
				}else{
					begin = 0;
				}
				if (a <= stop - 1000){
					end = a + 1000;
				}else{
					end = stop;
				}
				//above is new
				CloseableIterator<SAMRecord> iter = reader.query(chr, begin, end, false);
				double value = 0;
				while(iter.hasNext()){
					SAMRecord record = iter.next();
					if(!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()) {
						int tagstart = -1;
						int tagstop = -1;
					
						if(record.getInferredInsertSize() > 0) {
							tagstart = record.getAlignmentStart();
							
							tagstop = record.getAlignmentStart() + record.getInferredInsertSize() - 1;
						}else if (record.getInferredInsertSize() < 0 ) {
							tagstart = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
							tagstop = record.getAlignmentEnd();
						}
						if (tagstart == a || tagstop == a || (tagstart < a && tagstop > a)){
							value+=1;
						}
					}
					
				}
				iter.close();
				PileupNode2 node = new PileupNode2(a,value,chr);
				temp.add(node);
				mu.increment(value);
				std.increment(value);
			}
			PileupToBedGraph convert = new PileupToBedGraph(temp,1);
			bedgraph.addAll(convert.getBedGraph());
			
		}
		reader.close();
		mean = mu.getResult();
		stddev = std.getResult();
	}
	
}
