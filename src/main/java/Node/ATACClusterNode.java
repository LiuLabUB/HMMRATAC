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

public class ATACClusterNode extends ATACMatrixNode {
	
	
	private int _cluster;
	
	public ATACClusterNode(String chrom, int pos, double enrich1,
						   double enrich2, double enrich3, int index, int cluster) {
		super(chrom, pos, enrich1, enrich2, enrich3, index);
		_cluster = cluster;
	}
	
	public void setCluster(int cluster) {
		_cluster = cluster;
	}
	
	public int getCluster() {
		return _cluster;
	}
}
