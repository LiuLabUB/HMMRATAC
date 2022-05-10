package JAHMMTest;

import RobustHMM.ObservationReader;
import RobustHMM.TrainingReader;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TestDriver {
	
	
	private static String obs = null;
	private static File hmm = null;
	private static boolean supervise = false;
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {
			
			switch (args[i].charAt((1))) {
				
				case 'o':
					obs = (args[i + 1]);
					i++;
					break;
				case 't':
					String train = (args[i + 1]).toLowerCase();
					if (train.contains("t")) {
						supervise = true;
					}
					i++;
					break;
				case 'm':
					hmm = new File(args[i + 1]);
					i++;
					break;
				case 'h':
					printUsage();
					System.exit(0);
			}
		}
		if (obs == null || hmm == null) {
			printUsage();
			System.exit(0);
		}
		long start = System.nanoTime();
		List<List<ObservationVector>> newList = new ArrayList<>();
		if (!supervise) {
			ObservationReader obsReader = new ObservationReader(obs);
			List<?> obsList = obsReader.getObs();
			
			//split the data into 1000bp chunks
			int a;
			int i;
			int halfSize = obsList.size() / 2;
			for (i = 0; i < obsList.size() - 1; i += halfSize) {
				List<ObservationVector> temp = new ArrayList<>();
				for (a = i; a < i + halfSize - 1; a++) {
					
					ObservationVector o = (ObservationVector) obsList.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		} else {
			TrainingReader reader = new TrainingReader(obs);
			for (List<?> objects : reader.getObs()) {
				List<ObservationVector> temp = new ArrayList<>();
				for (Object object : objects) {
					temp.add((ObservationVector) object);
				}
				newList.add(temp);
			}
		}
		Hmm<?> h = null;
		if (hmm != null) {
			h = HmmBinaryReader.read(new FileInputStream(hmm));
		}
		
		Hmm<ObservationVector> firstHmm = (Hmm<ObservationVector>) h;
		checkModel(firstHmm);
		BaumWelchScaledLearner sbw = new BaumWelchScaledLearner();
		System.out.println("Initial Model\n" + firstHmm);
		KullbackLeiblerDistanceCalculator calc = new KullbackLeiblerDistanceCalculator();
		Hmm<ObservationVector> scaled = null;
		double distance;
		int iter = 0;
		while (iter < 100) {
			
			scaled = sbw.iterate(firstHmm, newList);
			checkModel(scaled);
			System.out.println("Model after iteration: " + iter + "\n" + scaled.toString());
			distance = (calc.distance(scaled, firstHmm) + calc.distance(firstHmm, scaled) / 2);
			System.out.println(distance);
			
			if (converged(scaled, firstHmm)) {
				break;
			}
			iter += 1;
			firstHmm = scaled;
		}
		
		File output = new File("model1.model");
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter.write(out, scaled);
		long end = System.nanoTime();
		long length = end - start;
		long seconds = length / 1000000000;
		System.out.println("****Time For Baum Welch =\t" + seconds + " seconds****");
	}
	
	public static void checkModel(Hmm<ObservationVector> hmm) {
		for (int i = 0; i < hmm.nbStates(); i++) {
			OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
			double[][] cov = pdf.covariance();
			FitRobust fitter = new FitRobust(cov);
			double[][] temp = fitter.getCovariance();
			OpdfMultiGaussian t = new OpdfMultiGaussian(pdf.mean(), temp);
			hmm.setOpdf(i, t);
		}
	}
	
	private static boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2) {
		int counter = 0;
		for (int i = 0; i < h1.nbStates(); i++) {
			OpdfMultiGaussian pdf1 = (OpdfMultiGaussian) h1.getOpdf(i);
			OpdfMultiGaussian pdf2 = (OpdfMultiGaussian) h2.getOpdf(i);
			double[] value1 = pdf1.mean();
			double[] value2 = pdf2.mean();
			for (int a = 0; a < value1.length; a++) {
				if (Math.abs(value1[a] - value2[a]) > 0.001) {
					counter++;
				}
			}
		}
		return counter == 0;
	}
	
	private static void printUsage() {
		System.out.println("Usage: java -jar ViterbiDriver.jar");
		System.out.println("Required Paramters:");
		System.out.println("-o <File> Observation File in specific format");
		System.out.println("Optional Paramters:");
		System.out.println("-t <boolean> Whether or not to perform supervised training");
		System.out.println("-m <File> File containing binary representation of HMM.");
	}
}
