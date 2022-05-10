package JAHMMTest;

import RobustHMM.ObservationReader;
import RobustHMM.TrainingReader;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import org.apache.commons.math3.stat.correlation.StorelessCovariance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class BasicHMM_EM {
	private static String obs = null;
	private static File model = null;
	private static boolean supervise = false;
	
	
	private final List<List<ObservationVector>> oseq;
	private Hmm<ObservationVector> hmm;
	
	public BasicHMM_EM(List<List<ObservationVector>> list, Hmm<ObservationVector> hmm) {
		this.oseq = list;
		this.hmm = hmm;
	}
	
	public Hmm<ObservationVector> getHMM() {
		return hmm;
	}
	
	public void iterate() {
		int dim = oseq.get(0).get(0).dimension();
		double[] initial = new double[hmm.nbStates()];
		double[][] trans = new double[hmm.nbStates()][hmm.nbStates()];
		ArrayList<StorelessCovariance> cov = new ArrayList<>();
		ArrayList<ArrayList<Mean>> mean = new ArrayList<>();
		for (int i = 0; i < hmm.nbStates(); i++) {
			StorelessCovariance temp = new StorelessCovariance(dim);
			cov.add(temp);
			ArrayList<Mean> t = new ArrayList<>();
			for (int a = 0; a < dim; a++) {
				Mean m = new Mean();
				t.add(m);
			}
			mean.add(t);
		}
		
		//Use current parameters to run viterbi
		for (List<ObservationVector> seq : oseq) {
			ViterbiCalculator vit = new ViterbiCalculator(seq, hmm);
			int[] states = vit.stateSequence();
			updateParameters(seq, states, initial, trans, cov, mean);
		}
		//Pi computation
		double iSum = Arrays.stream(initial).sum();
		for (int i = 0; i < initial.length; i++) {
			initial[i] /= iSum;
		}
		
		//Trans computation
		for (int i = 0; i < trans.length; i++) {
			int sum = 0;
			for (int a = 0; a < trans[i].length; a++) {
				sum += trans[i][a];
			}
			for (int a = 0; a < trans[i].length; a++) {
				trans[i][a] /= sum;
			}
		}
		List<OpdfMultiGaussian> opdf = new ArrayList<>();
		for (int i = 0; i < hmm.nbStates(); i++) {
			double[][] c = cov.get(i).getData();
			double[] m = new double[dim];
			for (int a = 0; a < mean.get(i).size(); a++) {
				m[a] = mean.get(i).get(a).getResult();
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(m, c);
			opdf.add(pdf);
		}
		hmm = new Hmm<>(initial, trans, opdf);
	}
	
	public void learn() {
		int maxIter = 100;
		for (int i = 0; i < maxIter; i++) {
			Hmm<ObservationVector> first = hmm;
			iterate();
			Hmm<ObservationVector> second = hmm;
			if (converged(first, second)) {
				break;
			}
		}
	}
	
	private void updateParameters(List<ObservationVector> seq, int[] states, double[] i,
								  double[][] t, ArrayList<StorelessCovariance> cov,
								  ArrayList<ArrayList<Mean>> mean) {
		i[states[0]]++;
		for (int a = 1; a < states.length; a++) {
			t[a - 1][a]++;
			
			cov.get(a).increment(seq.get(a).values());
			for (int x = 0; x < seq.get(a).dimension(); x++) {
				mean.get(a).get(x).increment(seq.get(a).value(x));
			}
		}
	}
	
	private boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2) {
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
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {
			
			switch (args[i].charAt((1))) {
				
				case 'o':
					obs = (args[i + 1]);
					i++;
					break;
				case 't':
					String train = (args[i + 1]);
					train = train.toLowerCase();
					if (train.contains("t")) {
						supervise = true;
					}
					i++;
					break;
				case 'm':
					model = new File(args[i + 1]);
					i++;
					break;
				case 'h':
					printUsage();
					System.exit(0);
			}
		}
		if (obs == null || model == null) {
			printUsage();
			System.exit(0);
		}
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
		Hmm<ObservationVector> h = null;
		if (model != null) {
			h = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(model));
		}
		BasicHMM_EM em = new BasicHMM_EM(newList, h);
		em.learn();
		Hmm<ObservationVector> outHMM = em.getHMM();
		File output = new File("model1.model");
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter.write(out, outHMM);
	}
	
	private static void printUsage() {
		System.out.println("Usage: java -jar ViterbiDriver.jar");
		System.out.println("Required Parameters:");
		System.out.println("-o <File> Observation File in specific format");
		System.out.println("Optional Parameters:");
		System.out.println("-t <boolean> Whether or not to perform supervised training");
		System.out.println("-m <File> File containing binary representation of HMM.");
	}
}
