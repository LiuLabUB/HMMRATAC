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
public class ORFNode {

	
		private String Chrom = "";
		private int Start = 0;
		private int Stop = 0;
		private int Direction = 0;
		private String ID = "";
		private double FPKM = 0.0;
		
		public ORFNode(){
			
		}
		
		public ORFNode(String chr,int start,int stop,int dir,String name,double fpkm)
		{
			Chrom = chr;
			Start = start;
			Stop = stop;
			Direction = dir;
			ID = name;
			FPKM = fpkm;
		}
		
		public double getFPKM(){
			return FPKM;
		}
		public void setFPKM(double fpkm){
			FPKM = fpkm;
		}
		
		public String getChrom(){
			return Chrom;
		}
		public void setChrom(String chr){
			Chrom = chr;
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
		
		public int getDir(){
			return Direction;
		}
		public void setDir(int dir){
			Direction = dir;
		}
		
		public String getID(){
			return ID;
		}
		public void setID(String name){
			ID = name;
		}

}
