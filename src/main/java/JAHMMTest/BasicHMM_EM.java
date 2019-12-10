package JAHMMTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.correlation.StorelessCovariance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import RobustHMM.ObservationReader;
import RobustHMM.TrainingReader;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;




public class BasicHMM_EM {
	private static String obs = null;
	private static File model = null;
	private static String train = null;
	private static boolean supervise = false;
	
	
	private List<List<ObservationVector>> oseq;
	private Hmm<ObservationVector> hmm;
	
	private int maxIter = 100;
	
	public BasicHMM_EM(List<List<ObservationVector>> list,Hmm<ObservationVector> h){
		oseq = list;
		hmm = h;
	}
	public Hmm<ObservationVector> getHMM(){return hmm;}
	public void iterate(){
		int dim = oseq.get(0).get(0).dimension();
		double[] initial = new double[hmm.nbStates()];
		double[][] trans = new double[hmm.nbStates()][hmm.nbStates()];
		ArrayList<StorelessCovariance> cov = new ArrayList<StorelessCovariance>();
		ArrayList<ArrayList<Mean>> mean = new ArrayList<ArrayList<Mean>>();
		for (int i = 0;i < hmm.nbStates();i++){
			StorelessCovariance temp = new StorelessCovariance(dim);
			cov.add(temp);
			ArrayList<Mean> t = new ArrayList<Mean>();
			for (int a = 0;a < dim;a++){
				Mean m = new Mean();
				t.add(m);
			}
			mean.add(t);
		}
		
		//Use current paramters to run viterbi
		for (List<ObservationVector> seq : oseq){
			ViterbiCalculator vit = new ViterbiCalculator(seq,hmm);
			int[] states = vit.stateSequence();
			updateParameters(seq,states,initial,trans,cov,mean);
		}
		//Pi computation
		int iSum=0;
		for (int i = 0;i < initial.length;i++){
			iSum += initial[i];
		}
		for (int i = 0;i < initial.length;i++){
			initial[i] /= iSum;
		}
		
		//Trans computation
		for (int i = 0; i < trans.length;i++){
			int sum = 0;
			for (int a = 0;a < trans[i].length;a++){
				sum += trans[i][a];
			}
			for (int a = 0;a < trans[i].length;a++){
				trans[i][a] /= sum;
			}
		}
		List<OpdfMultiGaussian> opdf = new ArrayList<OpdfMultiGaussian>();
		for (int i = 0; i < hmm.nbStates();i++){
			double[][] c = cov.get(i).getData();
			double[] m = new double[dim];
			for (int a = 0;a < mean.get(i).size();a++){
				m[a] = mean.get(i).get(a).getResult();
			}
			OpdfMultiGaussian pdf = new OpdfMultiGaussian(m,c);
			opdf.add(pdf);
		}
		hmm = new Hmm<ObservationVector>(initial,trans,opdf);
		
	}
	public void learn(){
		for (int i = 0; i < maxIter;i++){
			Hmm<ObservationVector> first = hmm;
			iterate();
			Hmm<ObservationVector> second = hmm;
			if (converged(first,second)){
				break;
			}
		}
	}
	private void updateParameters(List<ObservationVector> seq,int[] states,double[] i,
			double[][] t,ArrayList<StorelessCovariance> cov,
			ArrayList<ArrayList<Mean>> mean){
		i[states[0]]+=1;
		for (int a = 1;a < states.length;a++){
			if (a > 0){
				t[a-1][a]+=1;
			}
			
			cov.get(a).increment(seq.get(a).values());
			for (int x = 0;x < seq.get(a).dimension();x++){
				mean.get(a).get(x).increment(seq.get(a).value(x));
			}
		}
		
	}
	private  boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2){
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
	
	public static void main(String[] args) throws FileNotFoundException, IOException{
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
				model = new File(args[i+1]);
				i++;
				break;
			case'h':
				printUsage();
				System.exit(0);
			}
		}
		if (obs == null || model == null){
			printUsage();
			System.exit(0);	
		}
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
		Hmm<ObservationVector> h = null;
		if (model != null){
			h = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(model));
			//System.out.println(h.toString());
			//System.exit(0);
		}
		BasicHMM_EM em = new BasicHMM_EM(newList, h);
		em.learn();
		Hmm<ObservationVector> outHMM = em.getHMM();
		File output = new File("model1.model");
		FileOutputStream out = new FileOutputStream(output);
		HmmBinaryWriter writer = new HmmBinaryWriter();
		writer.write(out, outHMM);
		
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
