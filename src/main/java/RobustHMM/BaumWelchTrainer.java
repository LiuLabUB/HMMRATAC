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
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;

public class BaumWelchTrainer {
	
	private static String train = null;
	private static File hmm = null;
	private static int k = 0;
	private static int numDist = 0;
	private static String output=null;
	
	public static void main(String[] args) throws IOException{
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			case 't':
				train = (args[i+1]);
				i++;
				break;
			case'm':
				hmm = new File(args[i+1]);
				i++;	
				break;
			case'k':
				k = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'n':
				numDist = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'o':
				output = args[i+1];
				i++;
				break;
			}
		}
		
		if (train == null || (hmm == null && k == 0)){
			printUsage();
			System.exit(1);
		}
		
		TrainingReader reader = new TrainingReader(train);
		List<List<?>> trainList = reader.getObs();
		
		Hmm<?> h = null;
		if (hmm != null){
			h = HmmBinaryReader.read(new FileInputStream(hmm));
		}
		
		List<?> obs = trainList.get(0);
		String method = reader.getMethod();
		RobustHMM model = new RobustHMM(obs, trainList, h, true, k, method, numDist);
		Hmm<?> newHmm = model.getHmm();
		
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, newHmm);
		System.out.println(newHmm.toString());
	}

	public static void printUsage(){
		System.out.println("Usage: java -jar BaumWelchTrainer.jar");
		System.out.println("Parameters:");
		System.out.println("-t <File> Training file");
		System.out.println("-m <File> Optional Binary Model File");
		System.out.println("-k <int> Integer number of states for Kmeans Learning. Optional");
		System.out.println("-n <int> Integer number of distributions for mixture model or number of discreet integers for Integer model");
		System.out.println("-o <File> output of final model");
	}
}
