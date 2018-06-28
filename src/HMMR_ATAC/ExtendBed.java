package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class ExtendBed {
	
	private ArrayList<TagNode> input;
	private int extSize;
	private ArrayList<TagNode> output;
	
	/**
	 * Constructor for creating ExtendBed object and performing extension
	 * @param i an ArrayList of TagNode representing the BED data to be extended
	 * @param ext an integer representing the upstream and downstream extension size
	 */
	public ExtendBed(ArrayList<TagNode> i,int ext){
		input = i;
		extSize = ext;
		output = new ArrayList<TagNode>();
		set();
	}
	/**
	 * Access the extended data
	 * @return an ArrayList of TagNode representing the extended data
	 */
	public ArrayList<TagNode> getResults(){
		return output;
	}
	/**
	 * Extend the data by the extension size
	 */
	private void set(){
		
		for (int i = 0; i < input.size();i++){
			String chr = input.get(i).getChrom();
			int start = input.get(i).getStart();
			if (start - extSize > 0){
				start = start - extSize;
			}
			int stop = input.get(i).getStop() + extSize;
			TagNode temp = new TagNode(chr,start,stop);
			output.add(temp);
		}
	}
	

}
