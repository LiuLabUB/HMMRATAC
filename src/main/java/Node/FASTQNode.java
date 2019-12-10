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
public class FASTQNode {

	
		private String Header = "";
		private StringBuilder Sequence;
		private char Plus = '\0';
		private String Quality = "";
		
		public FASTQNode(){
			
		}
		
		public FASTQNode(String header,StringBuilder seq,char plus,String quality)
		{
			Header = header;
			Sequence = seq;
			Plus = plus;
			Quality = quality;
		}
		
		public String getHeader(){
			return Header;
		}
		public void setHeader(String header){
			Header = header;
		}
		
		public StringBuilder getSeq(){
			return Sequence;
		}
		public void setSeq(StringBuilder seq){
			Sequence = seq;
		}
		
		public char getPlus(){
			return Plus;
		}
		public void setPlus(char plus){
			Plus = plus;
		}
		
		public String getQuality(){
			return Quality;
		}
		public void setQuality(String quality){
			Quality = quality;
		}
}
