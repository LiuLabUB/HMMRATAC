package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class MergeBed {
	
	private ArrayList<TagNode> input;
	
	private ArrayList<TagNode> output;
	
	public MergeBed(ArrayList<TagNode> i){
		input = i;
		output = new ArrayList<TagNode>();
		merge();
	}
	public ArrayList<TagNode> getResults(){return output;}
	private void merge(){
		
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
