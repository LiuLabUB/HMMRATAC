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

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import GenomeFileReaders.bedFileReader;
import Node.TagNode;

public class WigCorrelation {

	
	private static String wig1 = null;
	private static String wig2 = null;
	private static String bed = null;
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			case'a':
				wig1 = args[i+1];
				i++;
				break;
			case'b':
				wig2 = args[i+1];
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
		
		if (wig1 == null || wig2 == null || bed == null){
			printUsage();
			System.exit(1);
		}
		bedFileReader reader = new bedFileReader(bed);
		ArrayList<TagNode> bedData = reader.getData();
		reader = null;
		BBFileReader wigReader1 = new BBFileReader(wig1);
		BBFileReader wigReader2 = new BBFileReader(wig2);
		ArrayList<double[]> mat = new ArrayList<double[]>();
		for (int i = 0; i < bedData.size();i++){
			String chrom = bedData.get(i).getChrom();
			int start = bedData.get(i).getStart();
			int stop = bedData.get(i).getStop();
			BigWigIterator iter1 = wigReader1.getBigWigIterator(chrom, start, chrom, stop, false);
			HashMap<Integer,double[]> map = new HashMap<Integer,double[]>();
			for (int a = start;a < stop;a++){
				double[] tmp = new double[2];
				map.put(a, tmp);
			}
			while(iter1.hasNext()){
				WigItem item = iter1.next();
				int begin = item.getStartBase();
				int end = item.getEndBase();
				double value = item.getWigValue();
				for (int a = begin;a < end;a++){
					if (a >= start && a < stop){
						double[] tmp = map.get(a);
						tmp[0] = value;
						map.put(a, tmp);
					}
				}
				
			}
			BigWigIterator iter2 = wigReader2.getBigWigIterator(chrom, start, chrom, stop,false);
			while(iter2.hasNext()){
				WigItem item = iter2.next();
				int begin = item.getStartBase();
				int end = item.getEndBase();
				double value = item.getWigValue();
				for (int a = begin;a < end;a++){
					if (a >= start && a < stop){
						double[] tmp = map.get(a);
						tmp[1] = value;
						map.put(a, tmp);
					}
				}
				
			}
			for (int a = start;a < stop;a++){
				mat.add(map.get(a));
			}
			
		}
		double[] matrix1 = new double[mat.size()];
		double [] matrix2 = new double[mat.size()];
		for (int i = 0;i < mat.size();i++){
			double[] value = mat.get(i);
			matrix1[i] = value[0];
			matrix2[i] = value[1];
		}
		mat = null;
		PearsonsCorrelation cor = new PearsonsCorrelation();
		double c = cor.correlation(matrix1,matrix2);
		matrix1 = null;matrix2 = null;
		System.out.println(wig1+"\t"+wig2+"\t"+c);
		
	}

	public static void printUsage(){
		System.out.println("Usage: java -jar WigCorrelation.jar <options>");
		System.out.println("Required Parameters:");
		System.out.println("\t-a <WIG/BigWIG> Wig File One");
		System.out.println("\t-b <WIG/BigWIG> Wig File Two");
		System.out.println("\t-r <BED> Regions of interest in BED format");
		System.out.println("\t-h Print this help message and exit");
		System.out.println("****IMPORTANT NOTE: java 7 required!!!****");
	}
}
