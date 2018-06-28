package HMMR_ATAC;

import java.util.ArrayList;
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
					int index;
					if (results.get(0).getHit().getStart() < inTemp.get(i).getStart()){
						output.add(new TagNode(inTemp.get(i).getChrom(),results.get(0).getHit().getStop(),
								results.get(1).getHit().getStart()));
						index = 1;
					}
					else{
						output.add(new TagNode(inTemp.get(i).getChrom(),inTemp.get(i).getStart(),
								results.get(0).getHit().getStart()));
						index = 0;
					}
					for (int x = index;x < results.size()-1;x++){
						output.add(new TagNode(inTemp.get(i).getChrom(),results.get(x).getHit().getStop()
								,results.get(x+1).getHit().getStart()));
					}
					if (results.get(results.size()-1).getHit().getStop() < inTemp.get(i).getStop()){
						output.add(new TagNode(inTemp.get(i).getChrom(),results.get(results.size()-1).getHit().getStop(),
								inTemp.get(i).getStop()));
					}
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

