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
public class HistogramBinner {
	
	private double BinStart = 0.0;
	private double BinStop = 0.0;
	private int BinNumber = 0;
	
	public HistogramBinner(){
		
	}

	public HistogramBinner(double binstart,double binstop,int binnumber){
		BinStart = binstart;
		BinStop = binstop;
		BinNumber = binnumber;
	}
	
	public double getStart(){
		return BinStart;
	}
	public void setStart(double binstart){
		BinStart = binstart;
	}
	
	public double getStop(){
		return BinStop;
	}
	public void setStop(double binstop){
		BinStop = binstop;
	}
	
	public int getNumber(){
		return BinNumber;
	}
	public void setNumber(int binnumber){
		BinNumber = binnumber;
	}
}
