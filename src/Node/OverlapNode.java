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
public class OverlapNode {

	
	private TagNode _node;
	private boolean _hit;
	private boolean _isConsumed;
	/**
	 * Constructor for creating OverlapNode object
	 * @param n a TagNode representing the overlaped portion
	 * @param h a boolean determining if the TagNode had an overlap
	 * @param c a boolean determining if the TagNode was completely overlapped
	 */
	public OverlapNode(TagNode n,boolean h,boolean c){
		_node = n;
		_hit = h;
		_isConsumed = c;
	}
	/**
	 * Access whether the TagNode was completely overlapped
	 * @return a boolean determining if the TagNode was completely overlapped
	 */
	public boolean isConsumed(){return _isConsumed;}
	/**
	 * Set whether the TagNode was completely overlapped
	 * @param c a boolean determining if the TagNode was completely overlapped
	 */
	public void setConsumed(boolean c){_isConsumed = c;}
	/**
	 * Access whether the TagNode had any overlap
	 * @return a boolean determining if the TagNode had any overlap
	 */
	public boolean hasHit(){return _hit;}
	/**
	 * Access the overlapped portion
	 * @return a TagNode representing the overlapped portion
	 */
	public TagNode getHit(){return _node;}
	/**
	 * Set whether the TagNode had any overlap
	 * @param c a boolean determining if the TagNode had any overlap
	 */
	public void setHit(boolean h){_hit = h;}
	/**
	 * Set the overlapped portion
	 * @param n a TagNode representing the overlapped portion
	 */
	public void setOverlap(TagNode n){_node = n;}
}
