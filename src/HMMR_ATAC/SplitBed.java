package HMMR_ATAC;

import java.util.ArrayList;

import Node.TagNode;

public class SplitBed {
	
	private ArrayList<TagNode> input;
	private ArrayList<TagNode> output;
	
	public SplitBed(ArrayList<TagNode> i){
		input = i;
		output = new ArrayList<TagNode>();
		split();
	}
	public ArrayList<TagNode> getResult(){return output;}
	private void split(){
		for (int i = 0; i < input.size();i++){
			String chrom = input.get(i).getChrom();
			int start = input.get(i).getStart();
			int stop = input.get(i).getStop();
			for (int x = start;x < stop;x+=25000000){
				int end = x + 25000000;
				if (end > stop){
					end = stop;
				}
				TagNode temp = new TagNode(chrom,x,end);
				output.add(temp);
			}
		}
	}

}
