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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;

public class BinaryModelReader {

	
	private static File input = null;
	public static void main(String[] args) throws FileNotFoundException, IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			case'i':
				input = new File(args[i+1]);
				i++;
				break;
			}
		}
		if (input == null){
			printUsage();
			System.exit(1);
		}
		Hmm<?> h = HmmBinaryReader.read(new FileInputStream(input));
		System.out.println(h.toString());
		//new
		for (int i = 0; i < h.nbStates();i++){
			
			if (i < h.nbStates()-1){System.out.print(h.getPi(i)+",");}
			else{System.out.println(h.getPi(i));}
		}
		for (int i = 0;i < h.nbStates();i++){
			for (int a = 0;a < h.nbStates();a++){
				if (a < h.nbStates()-1){
				System.out.print(h.getAij(i, a)+",");}
				else{System.out.println(h.getAij(i, a));}
			}
		}
		//new
		for (int i = 0;i < h.nbStates();i++){
			OpdfMultiGaussian opdf = (OpdfMultiGaussian) h.getOpdf(i);
			double[] means = opdf.mean();
			double[][] cov = opdf.covariance();
			for (int a = 0;a < means.length;a++){
				if (a != means.length-1){
					System.out.print(means[a]+",");
				}
				else{System.out.println(means[a]);}
			}
			printCovariance(cov);
		}
		/*
		int numStates = h.nbStates();
		for (int i = 0;i < numStates;i++){
			for (int a = 0; a < numStates;a++){
				System.out.print(h.getAij(i, a)+",");
			}
			System.out.println();
		}
		*/
	}
	public static void printCovariance(double[][] cov){
		for (int i = 0;i < cov.length;i++){
			for (int a = 0;a < cov[i].length;a++){
				if(a < cov[i].length-1){
					System.out.print(cov[i][a]+",");
				}
				else{
					System.out.println(cov[i][a]);
				}
			}
		}
	}
	
	private static void printUsage(){
		System.out.println("Usage: java -jar BinaryModelReader.jar -i BinaryModelFile");
	}
}
