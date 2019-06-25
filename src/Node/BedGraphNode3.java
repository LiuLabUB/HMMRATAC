package Node;


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
public class BedGraphNode3 extends BEDGraphWithIndexNode{

	
	private int _cluster;
	public BedGraphNode3(String chr, int start, int stop, double score,
			int index,int cluster) {
		super(chr, start, stop, score, index);
		_cluster = cluster;
		
	}

	public BedGraphNode3(String chr, int start, int stop, int score,
			int index,int cluster) {
		super(chr, start, stop, score, index);
		_cluster = cluster;
	}
	public BedGraphNode3(){
		
	}
	public int getCluster(){return _cluster;}
	public void setCluster(int c){_cluster = c;}
	
	
}
