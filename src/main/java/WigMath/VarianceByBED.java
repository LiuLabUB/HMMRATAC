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
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import GenomeFileReaders.bedFileReader;
import Node.TagNode;

public class VarianceByBED {
	private static String wig1 = null;
	
	private static String bed = null;
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			case'a':
				wig1 = args[i+1];
				i++;
				break;
			
			case'r':
				bed = args[i+1];
				i++;
				break;
			case'h':
				printUsage();
				System.exit(1);
			}
		}
		
		if (wig1 == null || bed == null){
			printUsage();
			System.exit(1);
		}

		bedFileReader reader = new bedFileReader(bed);
		ArrayList<TagNode> bedData = reader.getData();
		reader = null;
		BBFileReader wigReader1 = new BBFileReader(wig1);
		for (int i = 0; i < bedData.size();i++){
			String chrom = bedData.get(i).getChrom();
			int start = bedData.get(i).getStart();
			int stop = bedData.get(i).getStop();
			Variance var = new Variance();
			BigWigIterator iter1 = wigReader1.getBigWigIterator(chrom, start, chrom, stop, false);
			while(iter1.hasNext()){
				WigItem item = iter1.next();
				int begin = item.getStartBase();
				int end = item.getEndBase();
				double value = item.getWigValue();
				for (int a = begin;a < end;a++){
					var.increment(value);
				}
			}
			System.out.println(chrom+"\t"+start+"\t"+stop+"\t"+var.getResult());
		}
	}
	public static void printUsage(){
		System.out.println("Usage: java -jar WigPuller.jar <options>");
		System.out.println("Required Parameters:");
		System.out.println("\t-a <WIG/BigWIG> Wig File One");
		System.out.println("\t-r <BED> Regions of interest in BED format");
		System.out.println("\t-h Print this help message and exit");
		System.out.println("****IMPORTANT NOTE: java 7 required!!!****");
	}
}
