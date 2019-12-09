package WigMath;
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
import java.io.IOException;
import java.util.ArrayList;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class PullWigAboveCutoff {

	private ArrayList<TagNode> g;
	
	private String wig;
	private double cut;
	private ArrayList<TagNode> result;
	/**
	 * Constructor
	 * @param w a String representing the bigwig file
	 * @param c a double representing the zscored cutoff for reporting
	 * @param G an ArrayList of TagNode representing the regions to limit the search to
	 * @throws IOException
	 */
	public PullWigAboveCutoff(String w,double c, ArrayList<TagNode> G) throws IOException{
		wig = w;
		cut = c;
		g = G;
		set();
	}
	/**
	 * Access the regions that are above the cutoff
	 * @return
	 */
	public ArrayList<TagNode> getResults(){return result;}
	/**
	 * Set the data
	 * @throws IOException
	 */
	private void set() throws IOException{
		result = new ArrayList<TagNode>();
		
		for (int i = 0;i < g.size();i++){
			MeanAndStd results = new MeanAndStd(wig,g.get(i));
			double mean = results.getMean();
			double std = results.getStd();
			results=null;
			BBFileReader reader = new BBFileReader(wig);
			String chrom = g.get(i).getChrom();
			int begin = g.get(i).getStart();
			int end = g.get(i).getStop();
			BigWigIterator iter = reader.getBigWigIterator(chrom,begin,chrom,end,false);
			while (iter.hasNext()){
				WigItem item = iter.next();
				item.getChromosome();
				int start = item.getStartBase();
				int stop = item.getEndBase();
				double value = item.getWigValue();
				value = (value-mean) / std;
				if (value >= cut){
					TagNode temp = new TagNode(chrom,start,stop,value);
					result.add(temp);
				}
			
			}
			reader=null;
		}
		
	}
	
}
