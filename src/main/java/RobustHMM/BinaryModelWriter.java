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
import java.io.FileOutputStream;
import java.io.IOException;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscrete;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;

public class BinaryModelWriter {

	/*
	 * This case will take in a file with a specific format and write a binary HMM for future use.
	 */
	
	private static String input = null;
	private static File output = null;
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			case'i':
				input = (args[i+1]);
				i++;
				break;
			case 'o':
				output = new File(args[i+1]);
				i++;
				break;
			case'h':
				printUsage();
				System.exit(0);
			case'e':
				printExtendedHelp();
				System.exit(0);
			}
		}
		if (input == null || output == null){
			printUsage();
			//printExtendedHelp();
			System.out.println("Program Terminated");
			System.exit(0);
		}
		
		ModelFileReader reader = new ModelFileReader(input);
		Hmm<?> hmm = reader.getHmm();
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, hmm);
		System.out.println(hmm.toString());
	}

	private static void printUsage(){
		System.out.println("This program will take a textual representation of a Hidden Markov Model and output a binary representaion suitable for downstream tools"+"\n");
		System.out.println("Usage: java -jar BinaryModelWriter.jar"+"\n");
		System.out.println("Required Paramters:");
		System.out.println("-i <File> Input file with model parameters");
		System.out.println("-o <File> Output binary file containing model paramters. Suitable for input into downstream HMM programs");
		System.out.println("\n"+"Optional Parameters:"+"\n"+"-h Print this help message and exit");
		System.out.println("-e Print extended help message including details on format of input file");
	}
	private static void printExtendedHelp(){
		printUsage();
		System.out.println("Format of Input File:"+"\n");
		System.out.println("First line: ONE of Integer or Real or Vector or Mixture");
		System.out.println("Second Line: Comma separated list of Initial Probabilities for each state");
		System.out.println("Next N lines: Comma separated list of transition probabilities. N equals number of states. Each of N lines represents one line of transition matrix");
		System.out.println("Next N Lines: Observation Distributions. For Integer: each line is a comma separated list of probabilities for each integer observation");
		System.out.println("\t"+"For Gaussian (real) observations: Each line contains the mean,Variance of the distribution for that states");
		System.out.println("\t"+"For Mixture (mix of monovariate gaussians) Each states emission is represented as three consecutive lines: first line is comma sepaarted list of means. second line is list of variances, third line is list of weights");
		System.out.println("\t"+"For MultiVariateGaussians: First line represents the means of the distribution for a state. Next lines represent the covariance matrix for that states.  That order is repeated for each state");
	}
}
