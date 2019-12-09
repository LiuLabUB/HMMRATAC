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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import Node.OverlapNode;
import Node.TagNode;

public class SubtractBed {
	
	private ArrayList<TagNode> input;
	private ArrayList<TagNode> exclude;
	private ArrayList<TagNode> output;
	
	/**
	 * Constructor for creating a SubtractBed object and subtracting the data
	 * @param i an ArrayList of TagNode representing the data to be subtracted from
	 * @param e an ArrayList of TagNode representing the data to be subtracted by
	 */
	public SubtractBed(ArrayList<TagNode> i, ArrayList<TagNode> e){
		input = i;
		exclude = e;
		output = new ArrayList<TagNode>();
		subtract();
	}
	/**
	 * Access the subtracted data
	 * @return an ArrayList of TagNode representing the subtracted data
	 */
	public ArrayList<TagNode> getResults(){return output;}
	/**
	 * Subtract the data
	 */
	private void subtract(){
		HashMap<String,ArrayList<TagNode>> in = toMap(input);
		HashMap<String,ArrayList<TagNode>> ex = toMap(exclude);
		for (String chr : in.keySet()){
			if (ex.containsKey(chr)){
				
			
			ArrayList<TagNode> inTemp = in.get(chr);
			ArrayList<TagNode> exTemp = ex.get(chr);
			//TagNode temp = inTemp.get(0);
			for (int i = 0; i < inTemp.size();i++){
				ArrayList<OverlapNode> results = new ArrayList<OverlapNode>();
				for (int a = 0;a < exTemp.size();a++){
					OverlapNode node = overlap(inTemp.get(i),exTemp.get(a));
					if (node.hasHit()){
						// If A is completely consumed by B, dont add as there wont be a subtraction.
						if(!node.isConsumed()){
							results.add(node);
						}
					}
				}
				//No hits, therefore add the whole A entry into output
				if(results.size() == 0){
					output.add(inTemp.get(i));
				}
				//Only one hit recorded. Only need to find one set of outputs
				else if (results.size() == 1){
					TagNode hit = results.get(0).getHit();
					/*
					 * +++++++++++
					 *    ----
					 * ===    ====
					 */
					if (hit.getStart() > inTemp.get(i).getStart() && hit.getStop() < inTemp.get(i).getStop()){
						output.add(new TagNode(hit.getChrom(),inTemp.get(i).getStart(),hit.getStart()));
						output.add(new TagNode(hit.getChrom(),hit.getStop(),inTemp.get(i).getStop()));
					}
					/*
					 * ++++++++++
					 * -------
					 *        ===
					 */
					else if(hit.getStart() == inTemp.get(i).getStart()){
						output.add(new TagNode(hit.getChrom(),hit.getStop(),inTemp.get(i).getStop()));
					}
					/*
					 * +++++++++++
					 *    --------
					 * ===   
					 */
					else if(hit.getStop() == inTemp.get(i).getStop()){
						output.add(new TagNode(hit.getChrom(),inTemp.get(i).getStart(),hit.getStop()));
					}
					/*
					 *     +++++++++
					 * ---------
					 *     =====    
					 */
					else if (hit.getStart()<inTemp.get(i).getStart()){
						output.add(new TagNode(hit.getChrom(),inTemp.get(i).getStart(),hit.getStop()));
					}
					/*
					 * ++++++++
					 *     -------
					 *     ====
					 */
					else if (hit.getStart()>inTemp.get(i).getStart()){
						output.add(new TagNode(hit.getChrom(),hit.getStart(),inTemp.get(i).getStop()));
					}
				}
				//More than one overlap
				else if (results.size() > 1){
					
					
					//Scan the hits to look for which bases in A survive. Then report contiguous intervals that survive
					ArrayList<TagNode> res = new ArrayList<TagNode>();
					for (int y = 0;y< results.size();y++){
						res.add(results.get(y).getHit());
					}
					Collections.sort(res,  TagNode.basepairComparator);
					//Added above. changed below from results.get(i).getHit() to res.get(i) as tag node
					
					int index;
					if (res.get(0).getStart() < inTemp.get(i).getStart()){
						output.add(new TagNode(inTemp.get(i).getChrom(),res.get(0).getStop(),
								res.get(1).getStart()));
						index = 1;
					}
					else{
						output.add(new TagNode(inTemp.get(i).getChrom(),inTemp.get(i).getStart(),
								res.get(0).getStart()));
						index = 0;
					}
					for (int x = index;x < res.size()-1;x++){
						output.add(new TagNode(inTemp.get(i).getChrom(),res.get(x).getStop()
								,res.get(x+1).getStart()));
					}
					if (res.get(res.size()-1).getStop() < inTemp.get(i).getStop()){
						output.add(new TagNode(inTemp.get(i).getChrom(),res.get(res.size()-1).getStop(),
								inTemp.get(i).getStop()));
					}
					
					/**
					//New Approach. Use an array 
					int[] keep = new int[inTemp.get(i).getLength()];
					//String chr = inTemp.get(i).getChrom();
					int bedStart = inTemp.get(i).getStart();
					//int bedStop = inTemp.get(i).getStop();
					for (int x = 0;x < results.size();x++){
						int start = results.get(x).getHit().getStart();
						int stop = results.get(x).getHit().getStop();
						for (int y = start;y < stop;y++){
							if ((y-bedStart) >= 0 && (y-bedStart) < keep.length){
								keep[y-bedStart]++;
							}
						}
					}
					
					ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
					for (int x = 0;x < keep.length;x++){
						if(keep[i] == 0){
							PileupNode2 pnode = new PileupNode2(i+bedStart,0,chr);
							pile.add(pnode);
						}
					}
					output.addAll((new PileupToBedGraph(pile,1).getBedGraph()));
					**/
				}
				
			}
			
		}
			else{
				output.addAll(in.get(chr));
			}
		}
		
	}
	/**
	 * Determine if two entries overlap each other
	 * @param node1 a TagNode representing one entry
	 * @param node2 a TagNode representing a seconfd entry
	 * @return a OverlapNode representing the overlap between the two TagNode
	 */
	public static OverlapNode overlap(TagNode node1,TagNode node2){
		if(node1.getChrom().equals(node2.getChrom())){
			int start1 = node1.getStart();
			int start2 = node2.getStart();
			int stop1 = node1.getStop();
			int stop2 = node2.getStop();

			/*
			 * ++++++++++++
			 *    -----
			 *    =====
			 */
			if (start1 <= start2 && stop1 >= stop2){
				return new OverlapNode(node2,true,false);
			}
			
			/*
			 *    ++++++
			 *  -----  
			 *    ===
			 */
			else if(start1 >= start2 && start1 <= stop2){
				return new OverlapNode(new TagNode(node1.getChrom(),node1.getStart(),node2.getStop()),true,false);
			}
			
			/*
			 * ++++++++++
			 *      --------
			 *      =====
			 */
			else if (stop1 >= start2 && stop1 <= stop2){
				return new OverlapNode(new TagNode(node1.getChrom(),node2.getStart(),node1.getStop()),true,false);
			}
			
			/*
			 *    +++++++
			 * --------------   
			 */
			else if (start1 >= start2 && stop1 <= stop2){
				return new OverlapNode(node1,true,true);
			}
			else{
				return new OverlapNode(null,false,false);
			}
		} else{
			return new OverlapNode(null,false,false);
		}
	}
	/**
	 * Split the data by chromosome for calculation efficiency
	 * @param i an ArrayList of TagNode to split
	 * @return a HashMap of String and ArrayList of TagNode where the key String is the chromosome and the value ArrayList is all TagNode on that chromosome
	 */
	private HashMap<String,ArrayList<TagNode>> toMap(ArrayList<TagNode> i){
		HashMap<String,ArrayList<TagNode>> map = new HashMap<String,ArrayList<TagNode>>();
		for (int x = 0;x < i.size();x++){
			String chr = i.get(x).getChrom();
			if (map.containsKey(chr)){
				ArrayList<TagNode> temp = map.get(chr);
				temp.add(i.get(x));
				map.put(chr, temp);
			}
			else{
				ArrayList<TagNode> temp = new ArrayList<TagNode>();
				temp.add(i.get(x));
				map.put(chr, temp);
			}
		}
		
		return map;
		
	}
}

