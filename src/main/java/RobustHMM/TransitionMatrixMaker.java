package RobustHMM;
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
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import Node.TagNode;

public class TransitionMatrixMaker {

	
	private static String input = null;
	private static int numStates = 0;
	public static void main(String[] args) throws FileNotFoundException {
		for (int i = 0; i < args.length;i++){
			switch(args[i].charAt(1)){
			
			case'i':
				input = args[i+1];
				i++;
				break;
			case'n':
				numStates = Integer.parseInt(args[i+1]);
				i++;
				break;
			}
		}
		if (input == null || numStates == 0){
			printUsage();
			System.exit(1);
		}
		double initial = (double) (1.0/(double)numStates);
		for (int i = 0;i < numStates;i++){
			if (i < numStates-1){
				System.out.print(initial+",");
			}
			else{
				System.out.println(initial+",");
			}
		}
		
		StateReader reader = new StateReader(input);
		ArrayList<TagNode> list = reader.getList();
		double[][] trans = new double[numStates][numStates];
		for (int i = 0; i < list.size()-1;i++){
			int cluster = list.get(i).getScore(); ;
			int cluster2 = list.get(i+1).getScore();
			trans[cluster][cluster2]++;
			
		}
		for (int i = 0; i < trans.length;i++){
			
			int elementCounter=0;
			for (int x = 0;x < trans[i].length;x++){
				elementCounter += trans[i][x];
			}
			for (int y = 0; y < trans[i].length;y++){
				trans[i][y] /= elementCounter;
			}
		}

		for (int i = 0;i < trans.length;i++){
			for (int a = 0;a < trans[i].length;a++){
				if (a < trans[i].length-1){
					System.out.print(trans[i][a]+",");
				}
				else{
					System.out.println(trans[i][a]);
				}
			}
		}
	}
	
	public static void printUsage(){
		System.out.println("Usage: java -jar TransitionMatrixMaker.jar");
		System.out.println("-i <File> File containing state assignments. Format = chr\tstart\tstop\tstate");
		System.out.println("-n <int> Number of states");
	}
}

class StateReader{
	
	private ArrayList<TagNode> list;
	
	public StateReader(String input) throws FileNotFoundException{
		list = new ArrayList<TagNode>();
		TagNode temp = null;
		Scanner inFile = new Scanner((Readable) new FileReader(input));
		while (inFile.hasNext()){
			String line = inFile.nextLine();
			String[] feat = line.split("\\s+");
			String chr = feat[0];
			int start = Integer.parseInt(feat[1]);
			int stop = Integer.parseInt(feat[2]);
			int state = Integer.parseInt(feat[3]);
			for (int i = start;i < stop;i++){
				temp = new TagNode(chr,i,i+1,state);
				list.add(temp);
			}
			
		}
	}
	public ArrayList<TagNode> getList(){return list;}
}
