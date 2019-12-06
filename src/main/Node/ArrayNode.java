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

import java.util.Vector;


public class ArrayNode {

	private String Chrom = "";
	private int Start = 0;
	private int Stop = 0;
	private Vector<TagNode> Data = new Vector<TagNode>();
	private Vector<MyTagNode> List = new Vector<MyTagNode>();
		
	public ArrayNode(){
			
	}
	public ArrayNode(String chr,int start,int stop,Vector<TagNode> data,Vector<MyTagNode> list){
		Start = start;
		Stop = stop;
		Chrom = chr;
		Data = data;
		List = list;
	}
		public void setTemp(Vector<TagNode> data){
			Data = data;
		}
		
		public Vector<TagNode> getData(){
			return Data;
		}
		
		public void setList(Vector<MyTagNode> list){
			List = list;
		}
		public Vector<MyTagNode> getList(){
			return List;
		}
		
		
		public int getStart(){
			return Start;
		}
		public void setStart(int start){
			Start = start;
		}
		
		public int getStop(){
			return Stop;
		}
		public void setStop(int stop){
			Stop = stop;
		}
		
		public String getChrom(){
			return Chrom;
		}
		public void setChrom(String chr){
			Chrom = chr;
		}
	
	

}
