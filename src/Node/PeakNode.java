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
public class PeakNode {

	private int Summit = 0;
	private int Start = 0;
	private int Stop = 0;
	
	//public PeakNode(){
		
	//}
	
	public PeakNode(int summit,int start,int stop){
		Start = start;
		Summit = summit;
		Stop = stop;
	}
	
	public int getSummit(){
		return Summit;
	}
	public void setSummit(int summit){
		Summit = summit;
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
}
