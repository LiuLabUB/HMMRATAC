package stats;
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
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;

import GenomeFileReaders.bedFileReader;
import Node.TagNode;

public class BedOverlap {

	/*
	 * The idea is to feed one list into this program and call merge.
	 * This will return one list that is the overlapping regions of the first list.
	 * If you feed that merged list in again, it will return another merged list
	 */
	
	public BedOverlap(){
		
	}
	
	public ArrayList<TagNode> merge(ArrayList<TagNode> one,ArrayList<TagNode> two){
		ArrayList<TagNode> temp = new ArrayList<TagNode>();
		TagNode t = new TagNode();
		for (int i = 0;i < one.size();i++){
			TagNode node1 = one.get(i);
			for (int a = 0;a < two.size();a++){
				TagNode node2 = two.get(a);
				if (i != a){
					overlapNode result = intersect(node1,node2);
					if (result != null){
						t.setChrom(node1.getChrom());
						t.setStart(result.getStart());
						t.setStop(result.getStop());
						temp.add(t);
					}
				}
				
			}
		}
		
		return temp;
	}
	
	public overlapNode intersect(TagNode node1,TagNode node2){
		
		String chr = node1.getChrom();
		int start = node1.getStart();
		int stop = node1.getStop();
		String chr2 = node2.getChrom();
		int start2 = node2.getStart();
		int stop2 = node2.getStop();
		if (chr.equals(chr2)){
			if (start < start2 && stop > stop2){
				return new overlapNode(true,stop2-start2,start2,stop2);
			}
			else if (start > start2 && stop < stop2){
				return new overlapNode(true,stop-start,start,stop);
			}
			else if (start > start2 && stop > stop2 && start < stop2){
				return new overlapNode(true,stop2-start,start,stop2);
			}
			else if (start < start2 && stop < stop2 && stop < start2){
				return new overlapNode(true,stop-start2,start2,stop);
				
			}
			else{return null;}
		}
		else{return null;}
	}

	public class overlapNode{
		private boolean ans;
		private int overlap;
		private int overlapStart;
		private int overlapStop;
		public overlapNode(boolean a,int o,int start,int stop){
			ans = a;
			overlap = o;
			overlapStart = start;
			overlapStop = stop;
		}
		public void setStart(int start){overlapStart = start;}
		public int getStart(){return overlapStart;}
		public void setStop(int stop){overlapStop = stop;}
		public int getStop(){return overlapStop;}
		public void setAns(boolean a){
			ans = a;
		}
		public boolean getAns(){
			return ans;
		}
		public void setOverlap(int o){
			overlap = o;
		}
		public int getOverlap(){
			return overlap;
		}
	
	}
	private static String input;
	
	public static void main(String[] args){
		for (int i = 0 ; i < args.length;i++){
			switch(args[i].charAt(1)){
			case'i':
				input = args[i+1];
				i++;
				break;
			case'h':
				printUsage();
				System.exit(1);
			
			}
		}
		
		if (input == null ){
			printUsage();
			System.exit(1);
		}
		
		bedFileReader reader = new bedFileReader(input);
		ArrayList<TagNode> data = reader.getData();
		reader = null;
		
		BedOverlap over = new BedOverlap();
		for (int i = 0 ; i < 7; i++){
			data = over.merge(data, data);
		}
		
		for (int i = 0; i < data.size();i++){
			String chr = data.get(i).getChrom();
			int start = data.get(i).getStart();
			int stop = data.get(i).getStop();
			System.out.println(chr+"\t"+start+"\t"+stop);
		}
	}
	
	private static void printUsage(){
		System.out.println("Usage: java -jar BedOverlap.jar");
		System.out.println("Required Paramters:"+"\n"+"-i <BedFile> Input file");
		System.out.println("-h Print this help message and exit");
	}
}

