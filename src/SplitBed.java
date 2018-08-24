package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class SplitBed {
	
	private ArrayList<TagNode> input;
	private ArrayList<TagNode> output;
	private int window;
	
	/**
	 * Constructor to create new SplitBed object and split the data
	 * @param i an ArrayList of TagNode representing the data to be split
	 * @param w an integer representing the size of the windows to split the data into
	 */
	public SplitBed(ArrayList<TagNode> i,int w){
		input = i;
		window = w;
		output = new ArrayList<TagNode>();
		split();
	}
	/**
	 * Access the split data
	 * @return an ArrayList of TagNode representing the split data
	 */
	public ArrayList<TagNode> getResult(){return output;}
	/**
	 * Split the data by the window
	 */
	private void split(){
		for (int i = 0; i < input.size();i++){
			String chrom = input.get(i).getChrom();
			int start = input.get(i).getStart();
			int stop = input.get(i).getStop();
			for (int x = start;x < stop;x+=window){
				int end = x + window;
				if (end > stop){
					end = stop;
				}
				TagNode temp = new TagNode(chrom,x,end);
				output.add(temp);
			}
		}
	}

}
