package stats;
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

public class HG19Stats {
	
	private ArrayList<TagNode> stats = new ArrayList<TagNode>();
	public HG19Stats(){
		String[] s = new String[24];
		s[0] = "chr1	0	249250621";
		s[1] = "chr2	0	243199373";
		s[2] = "chr3	0	198022430";
		s[3] = "chr4	0	191154276";
		s[4] = "chr5	0	180915260";
		s[5] = "chr6	0	171115067";
		s[6] = "chr7	0	159138663";
		s[7] = "chrX	0	155270560";
		s[8] = "chr8	0	146364022";
		s[9] = "chr9	0	141213431";
		s[10] = "chr10	0	135534747";
		s[11] = "chr11	0	135006516";
		s[12] = "chr12	0	133851895";
		s[13] = "chr13	0	115169878";
		s[14] = "chr14	0	107349540";
		s[15] = "chr15	0	102531392";
		s[16] = "chr16	0	90354753";
		s[17] = "chr17	0	81195210";
		s[18] = "chr18	0	78077248";
		s[19] = "chr20	0	63025520";
		s[20] = "chrY	0	59373566";
		s[21] = "chr19	0	59128983";
		s[22] = "chr22	0	51304566";
		s[23] = "chr21	0	48129895";
		for (int i = 0;i < s.length;i++){
			String[] temp = s[i].split("\t");
			int start = Integer.parseInt(temp[1]);
			int stop = Integer.parseInt(temp[2]);
			TagNode node = new TagNode(temp[0],start,stop);
			stats.add(node);
			
		}
	}

	public ArrayList<TagNode> getStats(){return stats;}
}
