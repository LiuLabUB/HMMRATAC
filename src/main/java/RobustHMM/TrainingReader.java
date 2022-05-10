package RobustHMM;
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

import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class TrainingReader {
	
	private final String input;
	private List<List<?>> obs;
	private String method;
	
	public TrainingReader(String i) throws FileNotFoundException {
		input = i;
		read();
	}
	
	public List<List<?>> getObs() {
		return obs;
	}
	
	public String getMethod() {
		return method;
	}
	
	private void read() throws FileNotFoundException {
		Scanner inFile = new Scanner(new FileReader(input));
		obs = new ArrayList<>();
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			if (isInt(line)) {
				HashMap<Integer, ArrayList<ObservationInteger>> map = readInt(inFile);
				for (int index : map.keySet()) {
					List<?> tempList = map.get(index);
					obs.add(tempList);
				}
			} else if (isReal(line)) {
				HashMap<Integer, ArrayList<ObservationReal>> map = readReal(inFile);
				for (int index : map.keySet()) {
					List<?> tempList = map.get(index);
					obs.add(tempList);
				}
			} else if (isVector(line)) {
				HashMap<Integer, ArrayList<ObservationVector>> map = readVector(inFile);
				for (int index : map.keySet()) {
					List<?> tempList = map.get(index);
					obs.add(tempList);
				}
			}
		}
	}
	
	private HashMap<Integer, ArrayList<ObservationVector>> readVector(Scanner inFile) {
		HashMap<Integer, ArrayList<ObservationVector>> map = new HashMap<>();
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			int index = Integer.parseInt(temp[0]);
			double[] tempvalue = new double[temp.length - 1];
			for (int i = 1; i < temp.length; i++) {
				tempvalue[i - 1] = Double.parseDouble(temp[i]);
			}
			ObservationVector value = new ObservationVector(tempvalue);
			ArrayList<ObservationVector> list;
			list = map.containsKey(index) ? map.get(index) : new ArrayList<>();
			list.add(value);
			map.put(index, list);
		}
		return map;
	}
	
	private HashMap<Integer, ArrayList<ObservationReal>> readReal(Scanner inFile) {
		HashMap<Integer, ArrayList<ObservationReal>> map = new HashMap<Integer, ArrayList<ObservationReal>>();
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			int index = Integer.parseInt(temp[0]);
			double tempvalue = Double.parseDouble(temp[1]);
			ObservationReal value = new ObservationReal(tempvalue);
			ArrayList<ObservationReal> list;
			list = map.containsKey(index) ? map.get(index) : new ArrayList<>();
			list.add(value);
			map.put(index, list);
		}
		
		return map;
	}
	
	private HashMap<Integer, ArrayList<ObservationInteger>> readInt(Scanner inFile) {
		
		HashMap<Integer, ArrayList<ObservationInteger>> map = new HashMap<Integer, ArrayList<ObservationInteger>>();
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			int index = Integer.parseInt(temp[0]);
			int tempvalue = Integer.parseInt(temp[1]);
			ObservationInteger value = new ObservationInteger(tempvalue);
			ArrayList<ObservationInteger> list;
			list = map.containsKey(index) ? map.get(index) : new ArrayList<>();
			list.add(value);
			map.put(index, list);
		}
		return map;
		
	}
	
	private boolean isInt(String line) {
		if (line.toLowerCase().contains("integer") || line.toLowerCase().contains("int")) {
			method = "Integer";
			return true;
		} else {
			return false;
		}
	}
	
	private boolean isReal(String line) {
		if (line.toLowerCase().contains("real")) {
			method = "Real";
			return true;
		} else {
			return false;
		}
	}
	
	private boolean isVector(String line) {
		if (line.toLowerCase().contains("vector")) {
			method = "Vector";
			return true;
		} else {
			return false;
		}
	}
}
