package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class ExtendBed {
	
	private ArrayList<TagNode> input;
	private int extSize;
	private ArrayList<TagNode> output;
	
	public ExtendBed(ArrayList<TagNode> i,int ext){
		input = i;
		extSize = ext;
		output = new ArrayList<TagNode>();
		set();
	}
	public ArrayList<TagNode> getResults(){
		return output;
	}
	private void set(){
		
		for (int i = 0; i < input.size();i++){
			String chr = input.get(i).getChrom();
			int start = input.get(i).getStart();
			if (start - extSize > 0){
				start = start - extSize;
			}
			int stop = input.get(i).getStop() + extSize;
			TagNode temp = new TagNode(chr,start,stop);
			//System.out.println(chr+"\t"+start+"\t"+stop);
			output.add(temp);
		}
	}
	

}
