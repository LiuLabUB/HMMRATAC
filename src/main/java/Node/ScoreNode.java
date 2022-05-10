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

import java.util.ArrayList;

public class ScoreNode {
	private final double ave;
	private final double median;
	private final double max;
	private final ArrayList<Double> val;
	
	public ScoreNode(double ave, double med, double max, ArrayList<Double> val) {
		this.ave = ave;
		this.median = med;
		this.max = max;
		this.val = val;
	}
	
	public double getMean() {
		return ave;
	}
	
	public double getMax() {
		return max;
	}
	
	public double getMedian() {
		return median;
	}
	
	public ArrayList<Double> getvalues() {
		return val;
	}
}
