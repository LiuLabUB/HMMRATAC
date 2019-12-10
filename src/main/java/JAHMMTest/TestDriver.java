package JAHMMTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;
import RobustHMM.ModelFileReader;
import RobustHMM.ObservationReader;
import RobustHMM.TrainingReader;

public class TestDriver {

	
	private static String obs = null;
	private static File hmm = null;
	private static String train = null;
	private static boolean supervise = false;
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			case'o':
				obs = (args[i+1]);
				i++;
				break;
			case 't':
				train = (args[i+1]);
				train = train.toLowerCase();
				if(train.contains("t")){
					supervise = true;
				}
				i++;
				break;
			case'm':
				hmm = new File(args[i+1]);
				i++;
				break;
			case'h':
				printUsage();
				System.exit(0);
			}
		}
		if (obs == null || hmm == null){
			printUsage();
			System.exit(0);	
		}
		long start = System.nanoTime();
		List<List<ObservationVector>> newList = new ArrayList<List<ObservationVector>>();
		if (!supervise){
			ObservationReader obsReader = new ObservationReader(obs);
			List<?> obsList = obsReader.getObs();
			String method = obsReader.getMethod();
			
			//split the data into 1000bp chunks
			int a; int i;
			int halfSize = obsList.size()/2;
			//System.out.println(halfSize);
			for (i = 0;i < obsList.size()-1;i+=halfSize){
				//System.out.println("i is "+"\t"+i);
				List<ObservationVector> temp = new ArrayList<ObservationVector>();
				for (a = i;a < i+halfSize-1;a++){
					
					ObservationVector o = (ObservationVector) obsList.get(a);
					temp.add(o);
				}
				newList.add(temp);
			}
		}
		else{
			TrainingReader reader = new TrainingReader(obs);
			List<List<?>> temp = reader.getObs();
			for(int i = 0; i < temp.size();i++){
				List<ObservationVector> temp2 = new ArrayList<ObservationVector>();
				for(int a = 0;a < temp.get(i).size();a++){
					ObservationVector o = (ObservationVector) temp.get(i).get(a);
					temp2.add(o);
				}
				newList.add(temp2);
			}
		}
		Hmm<?> h = null;
		if (hmm != null){
			h = HmmBinaryReader.read(new FileInputStream(hmm));
			//System.out.println(h.toString());
			//System.exit(0);
		}
		
		Hmm<ObservationVector> firstHmm = (Hmm<ObservationVector>) h;
		firstHmm = checkModel(firstHmm);
		//BaumWelchLearner bw = new BaumWelchLearner();
		BaumWelchScaledLearner sbw = new BaumWelchScaledLearner();
		//BaumWelchScaledLearner sbw2 = new BaumWelchScaledLearner(2);
		System.out.println("Initial Model\n"+firstHmm.toString());
		//Hmm<ObservationVector> reg = bw.iterate(firstHmm, newList);
		//Hmm<ObservationVector> scaled = sbw.learn(firstHmm, newList);
		KullbackLeiblerDistanceCalculator calc = new KullbackLeiblerDistanceCalculator();
		Hmm<ObservationVector> scaled = null;
		double distance = 10;
		int iter = 0;
		while (iter < 100  ){
			
			scaled = sbw.iterate(firstHmm, newList);
			scaled = checkModel(scaled);
			System.out.println("Model after iteration: "+iter+"\n"+scaled.toString());
			distance = (calc.distance(scaled, firstHmm) + calc.distance(firstHmm, scaled) / 2);
			System.out.println(distance);
			
			if (converged(scaled,firstHmm)){
				break;
			}
			iter += 1;
			firstHmm = scaled;
		}
		
		//Hmm<ObservationVector> scaledWithConstant = sbw2.iterate(firstHmm, newList);
		
		//double[][] one = bw.getPoints();
		//ArrayList<Double> two = sbw.getPoints();
		//ArrayList<Double> three = sbw2.getPoints();
		
		//for (int b = 0;b < two.size();b++){
			
			//System.out.println(two.get(b)+"\t"+three.get(b));
		//}
		//System.out.println("Regular BW");
		//System.out.println(reg.toString());
		//System.out.println("Scaled BW");
		//System.out.println(scaled.toString());
		//System.out.println("Scaled with Constant");
		//System.out.println(scaledWithConstant.toString());
		
		File output = new File("model1.model");
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, scaled);
		long end = System.nanoTime();
		long length = end - start;
		long seconds = length / 1000000000;
		System.out.println("****Time For Baum Welch =\t"+seconds+" seconds****");
		/*
		output = new File("model2.model");
		FileOutputStream out2 = new FileOutputStream(output);
		HmmBinaryWriter writer2 = new HmmBinaryWriter();
		writer2.write(out2, scaledWithConstant);
		*/
	}
	public static Hmm<ObservationVector> checkModel(Hmm<ObservationVector> hmm){
		for (int i = 0;i < hmm.nbStates();i++){
			OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
			double[][] cov = pdf.covariance();
			FitRobust fitter = new FitRobust(cov);
			double[][] temp = fitter.getCovariance();
			OpdfMultiGaussian t = new OpdfMultiGaussian(pdf.mean(),temp);
			hmm.setOpdf(i, t);
		}
		return hmm;
	}
	private static boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2){
		int counter = 0;
		for (int i = 0;i < h1.nbStates();i++){
			OpdfMultiGaussian pdf1 = (OpdfMultiGaussian) h1.getOpdf(i);
			OpdfMultiGaussian pdf2 = (OpdfMultiGaussian) h2.getOpdf(i);
			double[] value1 = pdf1.mean();
			double[] value2 = pdf2.mean();
			for (int a = 0; a < value1.length;a++){
				//System.out.println((value1[a] - value2[a]));
				//System.out.println("Value1\t"+value1[a]+"\t"+"Value2\t"+value2[a]);
				if (Math.abs(value1[a] - value2[a]) > 0.001){
					
					counter += 1;
				}
			}
		}
		//System.out.println("counter\t"+counter);
		if (counter == 0){
			return true;
		}
		return false;
	}
	private static void printUsage(){
		System.out.println("Usage: java -jar ViterbiDriver.jar");
		System.out.println("Required Paramters:");
		System.out.println("-o <File> Observation File in specific format");
		System.out.println("Optional Paramters:");
		System.out.println("-t <boolean> Whether or not to perform supervised training");
		System.out.println("-m <File> File containing binary representation of HMM.");
	}

}
