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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;

public class ObservationReader {
	
	private String input;
	private List<?> obs;
	private String method;
	
	public ObservationReader(String i) throws FileNotFoundException{
		input = i;
		read();
	}

	public List<?> getObs(){
		return obs;
	}
	public String getMethod(){return method;}
	private void read() throws FileNotFoundException{
		Scanner inFile = new Scanner((Readable) new FileReader(input));
		
		while(inFile.hasNext()){
			String line = inFile.nextLine();
			//if (line.startsWith(">")){
				if (isInt(line)){
					method = "Integer";
					obs = readInt(inFile);
				}
				else if(isReal(line)){
					obs = readReal(inFile);
					method = "Real";
				}
				else if (isVector(line)){
					obs = readVector(inFile);
					method = "Vector";
				}
				else if(isMixture(line)){
					obs = readReal(inFile);
					method = "Mixture";
				}
			//}
			
		}
	}
	private List<?> readInt(Scanner inFile){
		List<ObservationInteger> obseq = new ArrayList<ObservationInteger>();
		while(inFile.hasNext()){
			ObservationInteger o = new ObservationInteger(inFile.nextInt());
			obseq.add(o);
		}
		return obseq;
	}
	private List<?> readReal(Scanner inFile){
		List<ObservationReal> obseq = new ArrayList<ObservationReal>();
		while(inFile.hasNext()){
			ObservationReal o = new ObservationReal(inFile.nextDouble());
			obseq.add(o);
		}
		return obseq;
	}
	private List<?> readVector(Scanner inFile){
		List<ObservationVector> obseq = new ArrayList<ObservationVector>();
		while(inFile.hasNext()){
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			double[] values = new double[temp.length];
			if (!line.contains("mask")){
				for (int i = 0; i < temp.length;i++){
					values[i] = Double.parseDouble(temp[i]);
				}
			}
			ObservationVector o = new ObservationVector(values);
			obseq.add(o);
		}
		return obseq;
	}
	
	private boolean isInt(String line){
		if (line.contains("Integer") || line.contains("integer") || line.contains("int") || line.contains("Int")){
			return true;
		}
		else{return false;}
	}
	private boolean isReal(String line){
		if (line.contains("Real") || line.contains("real")){
			return true;
		}
		else{return false;}
	}
	private boolean isVector(String line){
		if (line.contains("Vector") || line.contains("vector")){
			return true;
		}
		else{return false;}
	}
	private boolean isMixture(String line){
		if (line.contains("Mixture") || line.contains("mixture") || line.contains("mix") || line.contains("Mix")){
			return true;
		}
		else{return false;}
	}
	
}


