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
public class StateNode {

	private String chr;
	private int pos;
	private int state;
	
	public StateNode(String c, int p, int s){
		chr = c;
		pos = p;
		state = s;
	}
	public void setChr(String c){
		chr = c;
	}
	public String getChr(){
		return chr;
	}
	public void setPos(int p){
		pos = p;
	}
	public int getPos(){
		return pos;
	}
	public void setState(int s){
		state = s;
	}
	public int getState(){
		return state;
	}
	
}
