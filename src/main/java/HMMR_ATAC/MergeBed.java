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
