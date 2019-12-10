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
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;

public class ViterbiDriver {

	
	private static String obs = null;
	private static String train = null;
	private static File hmm = null;
	private static boolean welch = false;
	private static int k = 0;
	private static int numDist = 0;
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			case'o':
				obs = (args[i+1]);
				i++;
				break;
			case 't':
				train = (args[i+1]);
				i++;
				break;
			case'm':
				hmm = new File(args[i+1]);
				i++;
				break;
			case'w':
				String w = args[i+1].toLowerCase();
				if (w.equals("t") || w.equals("true")){welch = true;}
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
			case'h':
				printUsage();
				System.exit(0);
			}
		}
		if (obs == null || (train == null && hmm == null) ){
			printUsage();
			System.exit(0);
		}
		
		TrainingReader trainReader = null; List<List<?>> trainList = null;
		if (train != null){
			trainReader = new TrainingReader(train);
			trainList = trainReader.getObs();
		}
		ObservationReader obsReader = new ObservationReader(obs);
		List<?> obsList = obsReader.getObs();
		String method = obsReader.getMethod();
		trainReader=null;obsReader=null;
		Hmm<?> h = null;
		if (hmm != null){
			h = HmmBinaryReader.read(new FileInputStream(hmm));
			//System.out.println(h.toString());
			//System.exit(0);
		}
		
		RobustHMM HMM = new RobustHMM(obsList,trainList,h,welch,k,method,numDist);
		int[] states = HMM.getStates();
		Hmm<?> outHmm = HMM.getHmm();
		//System.out.println(outHmm.toString());
		for (int i = 0;i < states.length;i++){
			Observation o = (Observation) obsList.get(i);
			@SuppressWarnings("unchecked")
			Opdf<Observation> opdf = (Opdf<Observation>) outHmm.getOpdf(states[i]);
			double ep = Math.log10(opdf.probability(o));
			double tp;
			if (i == 0){
				tp = Math.log10(outHmm.getPi(states[i]));
			}
			else{
				tp = Math.log10(outHmm.getAij(states[i-1], states[i]));
			}
			double prob = ep + tp;
			System.out.println(o.toString()+"\t"+states[i]+"\t"+prob);
		}
		
		
	}

	private static void printUsage(){
		System.out.println("Usage: java -jar ViterbiDriver.jar");
		System.out.println("Required Paramters:");
		System.out.println("-o <File> Observation File in specific format");
		System.out.println("Optional Paramters:");
		System.out.println("-t <File> Training set file in specific format. Note if -w is true or -h is missing or -k is set, this is required");
		System.out.println("-m <File> File containing binary representation of HMM. If missing, -t is required");
		System.out.println("-w <T/F> Whether or not to use Baum Welch. If true, -t is required.");
		System.out.println("-k <int> Number of states used for KMeans training. If set, -t is required and -h is not required");
		System.out.println("-n <int> Number of Gaussian Distributions if using a mixture model or the number of discrete integers used for Integer model. Note: required if -k is set");
		System.out.println("-h Print this help message and exit.");
	}
}
