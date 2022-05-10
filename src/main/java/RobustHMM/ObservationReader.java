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
import java.util.List;
import java.util.Scanner;

public class ObservationReader {
	
	private final String input;
	private List<?> obs;
	private String method;
	
	public ObservationReader(String i) throws FileNotFoundException {
		input = i;
		read();
	}
	
	public List<?> getObs() {
		return obs;
	}
	
	public String getMethod() {
		return method;
	}
	
	private void read() throws FileNotFoundException {
		Scanner inFile = new Scanner(new FileReader(input));
		
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			if (isInt(line)) {
				method = "Integer";
				obs = readInt(inFile);
			} else if (isReal(line)) {
				obs = readReal(inFile);
				method = "Real";
			} else if (isVector(line)) {
				obs = readVector(inFile);
				method = "Vector";
			} else if (isMixture(line)) {
				obs = readReal(inFile);
				method = "Mixture";
			}
		}
	}
	
	private List<?> readInt(Scanner inFile) {
		List<ObservationInteger> obseq = new ArrayList<>();
		while (inFile.hasNext()) {
			ObservationInteger o = new ObservationInteger(inFile.nextInt());
			obseq.add(o);
		}
		return obseq;
	}
	
	private List<?> readReal(Scanner inFile) {
		List<ObservationReal> obseq = new ArrayList<>();
		while (inFile.hasNext()) {
			ObservationReal o = new ObservationReal(inFile.nextDouble());
			obseq.add(o);
		}
		return obseq;
	}
	
	private List<?> readVector(Scanner inFile) {
		List<ObservationVector> obseq = new ArrayList<>();
		while (inFile.hasNext()) {
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			double[] values = new double[temp.length];
			if (!line.contains("mask")) {
				for (int i = 0; i < temp.length; i++) {
					values[i] = Double.parseDouble(temp[i]);
				}
			}
			ObservationVector o = new ObservationVector(values);
			obseq.add(o);
		}
		return obseq;
	}
	
	private boolean isInt(String line) {
		return line.equalsIgnoreCase("integer") || line.equals("int");
	}
	
	private boolean isReal(String line) {
		return line.equalsIgnoreCase("real");
	}
	
	private boolean isVector(String line) {
		return line.equalsIgnoreCase("vector");
	}
	
	private boolean isMixture(String line) {
		return line.equalsIgnoreCase("mixture") || line.equals("mix");
	}
	
}


