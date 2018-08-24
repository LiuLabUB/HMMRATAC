package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class MergeBed {
	
	private ArrayList<TagNode> input;
	
	private ArrayList<TagNode> output;
	/**
	 * Constructor for creating new MergeBed object and merging data
	 * @param i an ArrayList of TagNode representing the data to be merged
	 */
	public MergeBed(ArrayList<TagNode> i){
		input = i;
		output = new ArrayList<TagNode>();
		merge();
	}
	/**
	 * Access the merged results
	 * @return an ArrayList of TagNode representing the merged data
	 */
	public ArrayList<TagNode> getResults(){return output;}
	/**
	 * Merge the data
	 */
	private void merge(){
		if (input.size()>0){
		TagNode first = input.get(0);
		for (int i = 1;i < input.size();i++){
			TagNode next = input.get(i);
			if (first.getStop() == next.getStart()){
				first.setStop(next.getStop());
			}
			else{
				output.add(first);
				first = next;
			}
		}
		output.add(first);
		}
	}

}
