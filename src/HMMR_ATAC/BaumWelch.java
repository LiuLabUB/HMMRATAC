package HMMR_ATAC;

import java.util.List;

import JAHMMTest.BaumWelchScaledLearner;
import JAHMMTest.FitRobust;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;

public class BaumWelch {
	
	private Hmm<?> h;
	private List<List<ObservationVector>> obs;
	
	public BaumWelch(Hmm<?> H, List<List<ObservationVector>> o){
		h = H;
		obs = o;
	}
	

	public Hmm<ObservationVector> build(){
		//System.out.println("Baum Welch started");
		Hmm<ObservationVector> firstHmm = (Hmm<ObservationVector>) h;
		firstHmm = checkModel(firstHmm);
		BaumWelchScaledLearner sbw = new BaumWelchScaledLearner();
		Hmm<ObservationVector> scaled = null;
		int iter = 0;
		while (iter < 100  ){
			//System.out.println("Baum welch iteration\t"+iter);
			scaled = sbw.iterate(firstHmm, obs);
			scaled = checkModel(scaled);
			if (converged(scaled,firstHmm)){
				break;
			}
			iter += 1;
			firstHmm = scaled;
		}
		//Set proportional initial probabilities
		for (int i = 0; i < scaled.nbStates();i++){
			double pi = (double)1/(double)scaled.nbStates();
			//System.out.println("Proportional Initial\t"+pi);
			scaled.setPi(i, 0.25);
		}
		return scaled;
	}
	
	private Hmm<ObservationVector> checkModel(Hmm<ObservationVector> hmm){
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
	
	private boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2){
		int counter = 0;
		for (int i = 0;i < h1.nbStates();i++){
			OpdfMultiGaussian pdf1 = (OpdfMultiGaussian) h1.getOpdf(i);
			OpdfMultiGaussian pdf2 = (OpdfMultiGaussian) h2.getOpdf(i);
			double[] value1 = pdf1.mean();
			double[] value2 = pdf2.mean();
			for (int a = 0; a < value1.length;a++){
				if (Math.abs(value1[a] - value2[a]) > 0.001){
					
					counter += 1;
				}
			}
		}
		
		if (counter == 0){
			return true;
		}
		return false;
	}
}
