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

import Node.OverlapNode;
import Node.TagNode;
import org.apache.commons.math3.ml.neuralnet.MapUtils;

import java.util.*;
import java.util.stream.Collectors;

public class SubtractBed {
	
	private final ArrayList<TagNode> input;
	private final ArrayList<TagNode> exclude;
	private final ArrayList<TagNode> output;
	
	/**
	 * Constructor for creating a SubtractBed object and subtracting the data
	 *
	 * @param input an ArrayList of TagNode representing the data to be subtracted from
	 * @param exclude an ArrayList of TagNode representing the data to be subtracted by
	 */
	public SubtractBed(ArrayList<TagNode> input, ArrayList<TagNode> exclude) {
		this.input = input;
		this.exclude = exclude;
		this.output = new ArrayList<>();
		subtract();
	}
	
	/**
	 * Access the subtracted data
	 *
	 * @return an ArrayList of TagNode representing the subtracted data
	 */
	public ArrayList<TagNode> getResults() {
		return output;
	}
	
	/**
	 * Subtract the data
	 */
	private void subtract() {
		Map<String, List<TagNode>> in = toMap(input);
		Map<String, List<TagNode>> ex = toMap(exclude);
		
		for (String chr : in.keySet()) {
			if (ex.containsKey(chr)) {
				
				List<TagNode> inTemp = in.get(chr);
				List<TagNode> exTemp = ex.get(chr);
				
				for (TagNode tagNode : inTemp) {
					ArrayList<OverlapNode> results = new ArrayList<>();
					for (TagNode node2 : exTemp) {
						OverlapNode node = overlap(tagNode, node2);
						if (node.hasHit()) {
							// If A is completely consumed by B, don't add as there won't be a subtraction.
							if (!node.isConsumed()) {
								results.add(node);
							}
						}
					}
					if (results.size() == 0) {
						//No hits, therefore add the whole A entry into output
						output.add(tagNode);
					} else if (results.size() == 1) {
						//Only one hit recorded. Only need to find one set of outputs
						TagNode hit = results.get(0).getHit();
						/*
						 * +++++++++++
						 *    ----
						 * ===    ====
						 */
						if (hit.getStart() > tagNode.getStart() && hit.getStop() < tagNode.getStop()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStart()));
							output.add(new TagNode(hit.getChrom(), hit.getStop(), tagNode.getStop()));
						}
						/*
						 * ++++++++++
						 * -------
						 *        ===
						 */
						else if (hit.getStart() == tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), hit.getStop(), tagNode.getStop()));
						}
						/*
						 * +++++++++++
						 *    --------
						 * ===
						 */
						else if (hit.getStop() == tagNode.getStop()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStop()));
						}
						/*
						 *     +++++++++
						 * ---------
						 *     =====
						 */
						else if (hit.getStart() < tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStop()));
						}
						/*
						 * ++++++++
						 *     -------
						 *     ====
						 */
						else if (hit.getStart() > tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), hit.getStart(), tagNode.getStop()));
						}
					} else {
						//More than one overlap
						
						//Scan the hits to look for which bases in A survive. Then report contiguous intervals that survive
						ArrayList<TagNode> res = new ArrayList<>();
						for (OverlapNode result : results) {
							res.add(result.getHit());
						}
						res.sort(TagNode.basepairComparator);
						
						int index;
						if (res.get(0).getStart() < tagNode.getStart()) {
							output.add(new TagNode(tagNode.getChrom(), res.get(0).getStop(),
									res.get(1).getStart()));
							index = 1;
						} else {
							output.add(new TagNode(tagNode.getChrom(), tagNode.getStart(),
									res.get(0).getStart()));
							index = 0;
						}
						for (int x = index; x < res.size() - 1; x++) {
							output.add(new TagNode(tagNode.getChrom(), res.get(x).getStop(),
									res.get(x + 1).getStart()));
						}
						if (res.get(res.size() - 1).getStop() < tagNode.getStop()) {
							output.add(new TagNode(tagNode.getChrom(), res.get(res.size() - 1).getStop(),
									tagNode.getStop()));
						}
						
						/*
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
			} else {
				output.addAll(in.get(chr));
			}
		}
	}
	
	/**
	 * Determine if two entries overlap each other
	 *
	 * @param node1 a TagNode representing one entry
	 * @param node2 a TagNode representing a second entry
	 * @return a OverlapNode representing the overlap between the two TagNode
	 */
	public static OverlapNode overlap(TagNode node1, TagNode node2) {
		if (node1.getChrom().equals(node2.getChrom())) {
			int start1 = node1.getStart();
			int start2 = node2.getStart();
			int stop1 = node1.getStop();
			int stop2 = node2.getStop();
			
			/*
			 * ++++++++++++
			 *    -----
			 *    =====
			 */
			if (start1 <= start2 && stop1 >= stop2) {
				return new OverlapNode(node2, true, false);
			}
			
			/*
			 *    ++++++
			 *  -----
			 *    ===
			 */
			else if (start1 >= start2 && start1 <= stop2) {
				return new OverlapNode(new TagNode(node1.getChrom(), node1.getStart(), node2.getStop()), true, false);
			}
			
			/*
			 * ++++++++++
			 *      --------
			 *      =====
			 */
			else if (stop1 >= start2 && stop1 <= stop2) {
				return new OverlapNode(new TagNode(node1.getChrom(), node2.getStart(), node1.getStop()), true, false);
			}
			
			/*
			 *    +++++++
			 * --------------
			 */
			else if (start1 >= start2 && stop1 <= stop2) {
				return new OverlapNode(node1, true, true);
			} else {
				return new OverlapNode(null, false, false);
			}
		} else {
			return new OverlapNode(null, false, false);
		}
	}
	
	/**
	 * Split the data by chromosome for calculation efficiency
	 *
	 * @param list an ArrayList of TagNode to split
	 * @return a Map of String and List of TagNode where the key String is the chromosome and the value List is all TagNodes on that chromosome
	 */
	private Map<String, List<TagNode>> toMap(ArrayList<TagNode> list) {
		return list.stream().collect(Collectors.groupingBy(TagNode::getChrom));
	}
}
