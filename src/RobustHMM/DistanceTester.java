package RobustHMM;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import JAHMMTest.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;

public class DistanceTester {

	private static File hmm = null;
	private static String train = null;
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws FileNotFoundException, IOException{
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			
			
			case'm':
				hmm = new File(args[i+1]);
				i++;
				break;
			case 't':
				train = (args[i+1]);
				i++;
				break;
			
			}
		}
		Hmm<?> h = null;
		h = HmmBinaryReader.read(new FileInputStream(hmm));
		TrainingReader trainReader = null; List<List<?>> trainList = null;
		if (train != null){
			trainReader = new TrainingReader(train);
			trainList = trainReader.getObs();
		}
		BaumWelchScaledLearner bw = new BaumWelchScaledLearner();
		Hmm<ObservationVector> first = (Hmm<ObservationVector>) h;
		List<List<ObservationVector>> newList = new ArrayList<List<ObservationVector>>();
		for (int i = 0;i < trainList.size();i++){
			ArrayList<ObservationVector> tempList = (ArrayList<ObservationVector>) trainList.get(i);
			newList.add(tempList);
		}
		int iter = 1;
		Hmm<ObservationVector> newHmm = first;
		newHmm = bw.iterate(newHmm, newList);
		
				
		//System.out.println("Number of iterations"+"\t"+iter);
		
		System.out.println("Number of iterations"+"\t"+bw.getNbIterations());

		KullbackLeiblerDistanceCalculator calc = new KullbackLeiblerDistanceCalculator();
		double dist1 = calc.distance(first, newHmm);
		double dist2 = calc.distance(newHmm, first);
		System.out.println("First Hmm");
		System.out.println(first.toString());
		System.out.println("Second HMM");
		System.out.println(newHmm.toString());
		System.out.println("distance1 = "+"\t"+dist1+"\tDistance2 ="+"\t"+dist2);
	}
}
