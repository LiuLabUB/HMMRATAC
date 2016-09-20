package RobustHMM;

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
