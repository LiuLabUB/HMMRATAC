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

import Node.TagNode;

import java.util.ArrayList;

public class SplitBed {
	
	private final ArrayList<TagNode> input;
	private final ArrayList<TagNode> output;
	private final int window;
	
	/**
	 * Constructor to create new SplitBed object and split the data
	 *
	 * @param input an ArrayList of TagNode representing the data to be split
	 * @param window an integer representing the size of the windows to split the data into
	 */
	public SplitBed(ArrayList<TagNode> input, int window) {
		this.input = input;
		this.window = window;
		this.output = new ArrayList<>();
		split();
	}
	
	/**
	 * Access the split data
	 *
	 * @return an ArrayList of TagNode representing the split data
	 */
	public ArrayList<TagNode> getResult() {
		return output;
	}
	
	/**
	 * Split the data by the window
	 */
	private void split() {
		for (TagNode tagNode : input) {
			String chrom = tagNode.getChrom();
			int start = tagNode.getStart();
			int stop = tagNode.getStop();
			for (int x = start; x < stop; x += window) {
				int end = Math.min(x + window, stop);
				output.add(new TagNode(chrom, x, end));
			}
		}
	}
}
