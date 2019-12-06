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
public class TemplateNode {
	
	private int _chrom;
	private int _pos;
	private char _dir;
	private int _count;
	
	public TemplateNode(int c, int p, char d, int ct){
		_chrom = c;
		_pos = p;
		_dir = d;
		_count = ct;
	}

	public void setChrom(int chr){_chrom = chr;}
	public int getChrom(){return _chrom;}
	public void setPos(int p){_pos = p;}
	public int getPos(){return _pos;}
	public void setDir(char d){_dir = d;}
	public char getDir(){return _dir;}
	public void setCount(int count){_count = count;}
	public int getCount(){return _count;}
	
}
