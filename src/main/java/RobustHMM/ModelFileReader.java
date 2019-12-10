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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixture;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;

public class ModelFileReader {
	
	private String input;
	private Hmm<?> hmm;
	
	private int numStates = 0;
	private double[] initial;
	private double[][] trans;
	//private List<? extends Opdf> opdf;
	
	public ModelFileReader(String file) throws FileNotFoundException{
		input = file;
		read();
	}
	public Hmm<?> getHmm(){
		return hmm;
	}
	
	private void read() throws FileNotFoundException{
		Scanner inFile = new Scanner((Readable) new FileReader(input));
		
		int counter = 0;
		while (inFile.hasNext()){
			String line = inFile.nextLine();
			if (counter == 0){
				if (isInt(line)){
					readInt(inFile);
				}
				else if (isReal(line)){
					readReal(inFile);
				}
				else if (isVector(line)){
					readVector(inFile);
				}
				else if (isMixture(line)){
					readMix(inFile);
				}
			}
		
			counter++;
		}
	}
	private void readMix(Scanner inFile){
		List<Opdf<ObservationReal>> opdf = new ArrayList<Opdf<ObservationReal>>();
		int counter = 0;
		int numVar = 0;
		ArrayList<String[]> t = new ArrayList<String[]>();
		while (inFile.hasNext()){
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			readInitialAndTrans(temp,counter);
			if (counter >= numStates+1){
				numVar = temp.length;
				t.add(temp);
			}
			
			counter++;
		}
		
		for (int i = 0; i < t.size();i+=3){
			double[] mean = new double[numVar];
			double[] var = new double[numVar];
			double[] prob = new double[numVar];
			for (int a = i;a < i+3;a++){
				String[] temp = t.get(a);
				for (int j = 0;j < temp.length;j++){
					double value = Double.parseDouble(temp[j]);
					if (a == i){
						mean[j] = value;
					}
					else if (a == i+1){
						var[j]=value;
					}
					else if (a == i+2){
						prob[j]=value;
					}
				}
			}
			
			OpdfGaussianMixture pdf = new OpdfGaussianMixture(mean,var,prob);
			opdf.add(pdf);
		}
		hmm = new Hmm<ObservationReal>(initial,trans,opdf);
	}
	private void readInt(Scanner inFile){
		List<Opdf<ObservationInteger>> opdf = new ArrayList<Opdf<ObservationInteger>>();
		int counter = 0;
		while(inFile.hasNext()){
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			readInitialAndTrans(temp,counter);
			
			OpdfInteger pdf = readIntPDF(temp,counter);
			if (pdf != null){
				opdf.add(pdf);
			}
			
			counter++;
		}
		hmm = new Hmm<ObservationInteger>(initial,trans,opdf);
	}
	
	private void readReal(Scanner inFile){
		List<Opdf<ObservationReal>> opdf = new ArrayList<Opdf<ObservationReal>>();
		int counter = 0;
		while(inFile.hasNext()){
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			readInitialAndTrans(temp,counter);
			OpdfGaussian pdf = readRealPDF(temp,counter);
			if (pdf != null){
				opdf.add(pdf);
			}
			
			
			counter++;
		}
		hmm = new Hmm<ObservationReal>(initial,trans,opdf);
	}
	private void readVector(Scanner inFile){
		int counter = 0;
		int numVar = 0;
		ArrayList<String[]> t = new ArrayList<String[]>();
		List<Opdf<ObservationVector>> opdf = new ArrayList<Opdf<ObservationVector>>();
		while(inFile.hasNext()){
			String line = inFile.nextLine();
			String[] temp = line.split(",");
			readInitialAndTrans(temp,counter);
			if (counter >= numStates+1){
				numVar = temp.length;
				t.add(temp);
			}
			
			
			counter++;
		}
		for (int i = 0;i< t.size();i+=(numVar+1)){
			double[] mean = new double[numVar];
			double[][] cov = new double[numVar][numVar];
			for (int a = i; a < (i+numVar+1);a++){
				if (a == i){
					for (int j = 0;j < t.get(a).length;j++){
						mean[j] = Double.parseDouble(t.get(a)[j]);
					}
				}
				else{
					for (int j = 0;j < t.get(a).length;j++){
						cov[a-(i+1)][j] = Double.parseDouble(t.get(a)[j]);
					}
				}
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(mean,cov);
			opdf.add(pdf);
		}
		hmm = new Hmm<ObservationVector>(initial,trans,opdf);
	}
	
	private boolean isInt(String line){
		if ( line.equals("Integer") || line.equals("integer") || line.equals("int")){
			return true;
		}
		else{
			return false;
		}
	}
	private boolean isReal(String line){
		if (line.equals("Real") || line.equals("real")){
			return true;
		}
		else{return false;}
	}
	private boolean isVector(String line){
		if (line.equals("Vector") || line.equals("vector")){
			return true;
		}
		else{return false;}
	}
	private boolean isMixture(String line){
		if (line.equals("Mixture") || line.equals("mixture") || line.equals("mix")){
			return true;
		}
		else{return false;}
	}
	private OpdfInteger readIntPDF(String[] temp,int counter){
		OpdfInteger pdf = null;
		if (counter > numStates){
			double[] probs = new double[temp.length];
			for (int i = 0 ; i < temp.length;i++){
				probs[i] = Double.parseDouble(temp[i]);
			}
			pdf = new OpdfInteger(probs);
			
		}
		return pdf;
	}
	private OpdfGaussian readRealPDF(String[] temp,int counter){
		OpdfGaussian pdf = null;
		if (counter > numStates){
			pdf = new OpdfGaussian(Double.parseDouble(temp[0]),Double.parseDouble(temp[1]));
		}
		
		return pdf;
	}
	
	
	private void readInitialAndTrans(String[] temp,int counter){
		if (counter == 0){
			numStates = temp.length;
			initial = new double[numStates];
			trans = new double[numStates][numStates];
			for (int i = 0;i < temp.length;i++){
				initial[i] = Double.parseDouble(temp[i]);
			}
		}
		if (counter > 0 && counter <= numStates){
			int i = counter - 1;
			for (int j = 0; j < numStates;j++){
				trans[i][j] = Double.parseDouble(temp[j]);
			}
		}
	}
	public static void main(String[] args) throws FileNotFoundException{
		ModelFileReader reader = new ModelFileReader("TestBEDFiles/TestMixtureModelFile");
		//reader.read();
		Hmm<?> hmm = reader.getHmm();
		//System.out.println(hmm.toString());
		Opdf<?> opdf = hmm.getOpdf(0);
		Observation obs = opdf.generate();
		System.out.println(obs.toString());
	}
}
