package HMMR_ATAC;

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

/*
 * Written by: Evan Tarbell 
 * evantarb@buffalo.edu
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;
import ATACFragments.FragPileupGen;
import FormatConverters.PileupToBedGraph;
import GEMM.HMMR_EM;
import GenomeFileReaders.GenomeFileReader;
import GenomeFileReaders.bedFileReader;
import Node.PileupNode2;
import Node.ScoreNode;
import Node.TagNode;
import RobustHMM.KMeansToHMM;
import RobustHMM.RobustHMM;

//import WigMath.FoldChange;//10_13_18
//import WigMath.PullWigAboveCutoff;//10_13_18

import WigMath.bedGraphMath;
import WigMath.pileup;

public class Main_HMMR_Driver {

	private static final boolean stopAfterModel = false;
	private static final boolean printHMMRTracks = false;

	/*
	 * Version number. Change as needed
	 */
	private static final String versionNum = "1.2.10";
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		
		ArgParser p = new ArgParser(args,versionNum);

		//Required inputs

		File bam = p.getBam();
		File index1 = p.getIndex();
		String genomeFile = p.getGenome();

		//Optional inputs

		//comma separated list of initial mean values for frag dist
		String means = p.getMeans();
		//comma separated list of initial standard deviations for frag dist
		String stddevs = p.getStd();
		//whether to perform fragment dist em
		boolean fragEM = p.getEM();
		//minimum mapping quality of reads to keep
		int minMapQ = p.getMinQ();
		//lower bound for fold change range for choosing training sites
		int lower = p.getLower();
		//upper bound for fold change range for choosing training sites
		int upper = p.getUpper();
		//zscored read coverage to exclude from viterbi decoding
		int zscore = p.getZscore();
		//output name
		String output = p.getOutput();
		// whether to print peaks
		boolean peaks = p.getPeaks();
		// whether to print bedgraph
		boolean bg = p.getBedgraph();
		String blacklist = p.getBlacklist();
		int minLength = p.getMinLength();
		String scoreSys = p.getScore();
		boolean BGScore = p.getBGScore();
		int k = p.getK();
		int trim = p.getTrim();
		String trainingRegions = p.getTrainingRegions();
		int vitWindow = p.getWindow();
		File modelFile = p.getModelFile();
		int maxTrain = p.getMaxTrain();
		boolean rmDup = p.getRemoveDuplicates();
		boolean printExclude = p.getPrintExclude();
		boolean printTrain = p.getPrintTrain();
		long randomTrainSeed = p.getRandomTrainSeed();
		double threshold = p.getThreshold();
		//For run time calculation
		long startTime = System.currentTimeMillis();
		
		//Declare output name
		if (output ==null){
			output = "NA";
		}
		PrintStream log = new PrintStream(output +".log");
		
		//Exit program if BAM file or Index file not given
		if (bam == null || index1 == null || genomeFile == null ){
			p.printUsage();
			System.exit(1);
		}
		
		/*
		 * Report version into log file
		 * 
		 */
		log.println("Version:"+"\t"+versionNum);
		
		/*
		 * Report all input arguments into log file
		 */
		log.println("Arguments Used:");
		for (int i = 0; i < args.length-1;i+=2){
			log.println(args[i]+"\t"+args[i+1]);
		}
		
		//Read in genome size stats
		GenomeFileReader gReader = new GenomeFileReader(genomeFile);
		ArrayList<TagNode> genomeStats = gReader.getMap();
		
		//Read in blacklisted if inputted
		ArrayList<TagNode> black = null;
		if (blacklist != null){
			black = new bedFileReader(blacklist).getData();
			
		}
		
		/*
		 * Set fragment length distribution parameters. 
		 * Use inputted values to set initial values, if provided.
		 * Else use defaults
		 */
		
		double[] fragMeans = new double[4];
		double[] fragStddevs = new double[4];
		double[] mode = new double[4];
		mode[1] = mode[2] = mode[3] = 2;
		mode[0]=0.5;
		
		if (means != null){
			String[] mu = means.split(",");
			for (int i = 0; i < mu.length;i++){
				fragMeans[i] = Double.parseDouble(mu[i]);
			}
		} else{
			fragMeans[0] = 50.0;
			fragMeans[1] = 200.0; 
			fragMeans[2] = 400.0;
			fragMeans[3] = 600.0;
		}
		if (stddevs != null){
			String[] std = stddevs.split(",");
			for(int i = 0;i < std.length;i++){
				fragStddevs[i] = Double.parseDouble(std[i]);
			}
		} else{
			Arrays.fill(fragStddevs, 20.0);
		}
		
		
		/*
		 * Pull the lengths from the read data. Fragment distribution uses fragments with length > 100 to train the 
		 nucleosome distributions. the short distribution is set to 50 and remains unchanged
		 Only occurs if EM training occurs, else use the default settings
		*/
		
		
		if (fragEM){
			pullLargeLengths puller = new pullLargeLengths(bam, index1, minMapQ, genomeStats,java.util.Arrays.copyOf(fragMeans, 4));
			double[] lengths = puller.getSampledLengths(10);
			double[] weights = puller.getWeights();

			//Perform EM training
			
			HMMR_EM em = new HMMR_EM(weights,java.util.Arrays.copyOfRange(fragMeans, 1,4),
					java.util.Arrays.copyOfRange(fragStddevs, 1, 4),lengths);
			em.learn();
			double[] tempMeans = em.getMeans();
			double[] tempLam = em.getLamda();
			for (int i = 0;i < tempMeans.length;i++){
				//This will update the parameters IFF they were updated. If they become NaN, leave as default
				if(!Double.isNaN(tempMeans[i]) && !Double.isNaN(tempLam[i])){
					fragMeans[i+1] = tempMeans[i];
					fragStddevs[i+1] = tempLam[i];
				}
			}
		}
		
		log.println("Fragment Expectation Maximum Done");
		for (int i = 0;i < fragMeans.length;i++){
			log.println("Mean\t"+fragMeans[i]+"\tStdDevs\t"+fragStddevs[i]);
		}
		
		
		/*
		 * Generate genome wide pileup of read coverage and simultaneously calculate mean and standard deviation
		 * to calculate fold change and zscore. Convert pileup to bedgraph for easier storage. Calculate the pileup
		 * in chunks for memory efficiency
		 * 
		 * NOTE: this method does not currently exist. Instead, we use an inputted big wig file to find fold change
		 * training regions and zscored masked regions
		 */
		
		
		
		/*
		 * Below is method to use BigWig input file
		 * Find regions with fold change within determined range to use as training sites.
		 * Find regions with zscore values above certain cutoff to exclude from viterbi.
		 */
		Hmm<ObservationVector> hmm=null;

		
		pileup pileupData = new pileup(new SplitBed(genomeStats, vitWindow).getResult(), 0, bam, index1, 0, rmDup);
		bedGraphMath fc = new bedGraphMath(pileupData.getBedGraph());
		
		//calculate the cpm scaling factor for input into FragPileupGen. Use cpmScale=1 for no scaling
		double cpmScale = pileupData.getCPMScale()/1000000;
		log.println("ScalingFactor\t"+cpmScale);
		
		double genomeMean = fc.getMean();
		double genomeStd = fc.getSTD();
		
		ArrayList<TagNode> train = new MergeBed(fc.getBetweenRanges(upper, lower)).getResults();

		
		ArrayList<TagNode> newTrain = new ArrayList<>();
		if (train.size() < maxTrain){
			maxTrain = train.size();
		}
		
		//Shuffle training list before choosing.
		Collections.shuffle(train, new Random(randomTrainSeed));
		for (int i = 0; i < maxTrain; i++){
			newTrain.add(train.get(i));
		}
		
		train = newTrain;
		train = new ExtendBed(train,5000).getResults();
		ArrayList<TagNode> exclude = new MergeBed(fc.getAboveZscore(zscore)).getResults();
		ArrayList<TagNode> addBack = exclude;

		if (blacklist != null){
			exclude.addAll(black);
			exclude = new MergeBed(exclude).getResults();
		}
			
		if(printExclude){
			PrintStream ex = new PrintStream(output +"_excluded.bed");
			for (TagNode tagNode : exclude) {
				ex.println(tagNode.toString() + "\t" + "exclude");
			}
			ex.close();
		}
		
		
		newTrain = new ArrayList<>();
		for (TagNode tagNode : train) {
			int counter = 0;
			for (TagNode node2 : exclude) {
				if (SubtractBed.overlap(tagNode, node2).hasHit()) {
					counter++;
				}
			}
			if (counter == 0) {
				newTrain.add(tagNode);

			}
		}
		
		train = newTrain;
		// Allows user to use training set for model generation
		if (trainingRegions != null){
			bedFileReader trainReader = new bedFileReader(trainingRegions);
			train = trainReader.getData();
		}

		log.println("Training Regions found and Zscore regions for exclusion found");
		
		
		/*
		 * Create the fragment pileup tracks using the training set and the fragment distribution parameters
		 */
		if (modelFile == null){
		
			
			if (printTrain) {
				PrintStream tr = new PrintStream(output + "_training.bed");
				for (TagNode tagNode : train) {
					tr.println(tagNode.toString() + "\t" + "training");
				}
				tr.close();
			}
			FragPileupGen gen = new FragPileupGen(bam, index1, train, mode, fragMeans, fragStddevs, minMapQ, rmDup,cpmScale);
			TrackHolder holder = new TrackHolder((gen.transformTracks(gen.scaleTracks(gen.getAverageTracks()))), trim);

			
			log.println("Training Fragment Pileup completed");
			
			/*
			 * Create the initial model using KMeans and then refine it using Baum-Welch
			 */
			
			KMeansToHMM kmeans = new KMeansToHMM(holder.getDataSet(), k,Integer.MAX_VALUE,true,true,true);
			log.println("Kmeans Model:\n"+kmeans.getHMM().toString());

			hmm = new BaumWelch(kmeans.getHMM(),holder.getBWObs(),1000).build();
		
		}
		/*
		 * Use input model if available
		 */
		
		if (modelFile != null){
			hmm = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(modelFile));
		}

		
		/*
		 * Identify peak state as the state with the highest short read signal.
		 * Identify flanking nucleosome state as the state with the second-highest mono read signal.
		 */
		
		int peak = -1;
		double max = 0.0;
		for (int i = 0; i < hmm.nbStates();i++){
			OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
			double sh = pdf.mean()[0];
			if (sh > max){
				peak = i;
				max = sh;
			}
		}
		max = 0.0;
		for (int i = 0; i < hmm.nbStates();i++){
			if (i != peak){
				OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
				double mono = pdf.mean()[1];
				if (mono > max){
					max = mono;
				}
			}
		}
		
		/*
		 * Output binary model file
		 */
		File outputModel = new File(output +".model");
		FileOutputStream outModel = new FileOutputStream(outputModel);
		HmmBinaryWriter.write(outModel, hmm);
		outModel.close();
		log.println("Model created and refined. See "+ output +".model");
		log.println("Model:\n"+ hmm);
		
		/*
		 * Stop program if only model is desired
		 */
		
		if(stopAfterModel){
			System.exit(0);
		}
		
		/*
		 * Split the genome file into smaller 25MB chunks 
		 * Can also split into whatever sized chunks the users prefers
		 * May be necessary to split into smaller chunks for machines with less memory
		 */
		
		ArrayList<TagNode> split = new SplitBed(genomeStats, vitWindow).getResult();
		
		
		/*
		 * Subtract excluded regions from the split genome for Viterbi
		 */
		
		ArrayList<TagNode> vitBed = new SubtractBed(split,exclude).getResults();

		log.println("Genome split and subtracted masked regions");
		
		/*
		 * Run viterbi on the whole genome
		 */
		PrintStream NFR = null;
		PrintStream MONO = null;
		PrintStream DI=null;
		PrintStream TRI=null;
		if(printHMMRTracks){
			 NFR = new PrintStream(output +"_nfr.bedgraph");
			 MONO = new PrintStream(output +"_mono.bedgraph");
			 DI = new PrintStream(output +"_di.bedgraph");
			 TRI = new PrintStream(output +"_tri.bedgraph");
		}
		ArrayList<TagNode> genomeAnnotation = new ArrayList<>();
		for (int i = 0;i < vitBed.size();i++){
			if (vitBed.get(i).getLength() >= 10){
				ArrayList<TagNode> tempBed = new ArrayList<>();
				tempBed.add(vitBed.get(i));
				FragPileupGen vGen = new FragPileupGen(bam, index1, tempBed, mode, fragMeans, fragStddevs, minMapQ, rmDup,cpmScale);
				TrackHolder vHolder = new TrackHolder(vGen.transformTracks(vGen.scaleTracks(vGen.getAverageTracks())), trim);
				
				if (printHMMRTracks){
					HMMRTracksToBedgraph tracks = new HMMRTracksToBedgraph(vHolder.getRawData(),vitBed.get(i),10);
					ArrayList<TagNode> nfr = tracks.getShort();
					ArrayList<TagNode> mono = tracks.getMono();
					ArrayList<TagNode> di = tracks.getDi();
					ArrayList<TagNode> tri = tracks.getTri();
					
					
					if (nfr != null) {
						for (TagNode tagNode : nfr) {
							NFR.println(tagNode.toString2());
						}
					}
					if (mono != null) {
						for (TagNode tagNode : mono) {
							MONO.println(tagNode.toString2());
						}
					}
					if (di != null) {
						for (TagNode tagNode : di) {
							DI.println(tagNode.toString2());
						}
					}
					if (tri != null) {
						for (TagNode tagNode : tri) {
							TRI.println(tagNode.toString2());
						}
					}
					
				}

				RobustHMM HMM = new RobustHMM(vHolder.getObs(),null,hmm,false,0,"Vector",0);
				int[] states = HMM.getStates();
				
				int start = vitBed.get(i).getStart();
				int remainder = vitBed.get(i).getLength() % 10;
				ArrayList<PileupNode2> pile = new ArrayList<>();
				int a;
				for (a = 0;a < states.length-1;a++){
					PileupNode2 pNode = new PileupNode2(start+(a*10),(double)states[a],vitBed.get(i).getChrom());
					pile.add(pNode);
				}
				PileupNode2 pNode = new PileupNode2(start+(((a)*10)-remainder),(double)states[a],vitBed.get(i).getChrom());
				pile.add(pNode);
				genomeAnnotation.addAll(new PileupToBedGraph(pile,10).getBedGraph());
				
				if(i%50 == 0 || i == vitBed.size()-1){
					log.println(i+" round viterbi done");
				}
			}
		}
		//out.close();
		if(printHMMRTracks){
			NFR.close();MONO.close();DI.close();TRI.close();
		}
		/*
		 * Report the final results as peaks, bedgraphs and summits, if desired
		 */
		PrintStream bedgraph=null;
		if (bg){
			 bedgraph = new PrintStream(output +".bedgraph");
		}
		PrintStream pks=null;
		PrintStream summits=null;
		if (peaks){
			 pks = new PrintStream(output +"_peaks.gappedPeak");
			 summits = new PrintStream(output +"_summits.bed");
		}
		HashMap<String,ArrayList<TagNode>> bdg = fc.getMappedBedgraph();

		HashMap<String,ArrayList<TagNode>> hmmrBdg = bedGraphMath.toMap(genomeAnnotation);
		int counter=1;
		for (String chr : hmmrBdg.keySet()){
			ArrayList<TagNode> hmmr = hmmrBdg.get(chr);
			ArrayList<TagNode> signal = bdg.get(chr);
			if (signal != null) {
				hmmr.sort(TagNode.basepairComparator);
				
				signal.sort(TagNode.basepairComparator);
				
				int index = 0;
				for (int i = 0; i < hmmr.size(); i++) {
					TagNode temp = hmmr.get(i);

					/*
					 * Execute the scoring commands if the state is a peak or if bedgraph scoring is on
					 */
					if ((int) temp.getScore2() == peak || BGScore) {
						boolean hasHadOverlap = false;
						ArrayList<TagNode> overlaps = new ArrayList<>();
						for (int a = index; a < signal.size(); a++) {
							if (SubtractBed.overlap(temp, signal.get(a))
									.hasHit()) {
								overlaps.add(signal.get(a));
								hasHadOverlap = true;
							} else {
								if (hasHadOverlap) {
									index = a;
									break;
								}
							}

						}
						ScoreNode scores = bedGraphMath.set(temp, overlaps);
						double MAX = scores.getMax();
						double MEAN = scores.getMean();
						double MEDIAN = scores.getMedian();
						double ZSCORE = (scores.getMean() - genomeMean) / genomeStd;
						double FOLDCHANGE = scores.getMean() / genomeMean;

						switch (scoreSys) {
							case "ave":
								temp.setScore3(Double.toString(MEAN));
								break;
							case "fc":
								temp.setScore3(Double.toString(FOLDCHANGE));
								break;
							case "zscore":
								temp.setScore3(Double.toString(ZSCORE));
								break;
							case "med":
								temp.setScore3(Double.toString(MEDIAN));
								break;
							case "all":
								String ANSWER = MAX + "_" + MEAN + "_" + MEDIAN + "_" + ZSCORE + "_" + FOLDCHANGE;
								temp.setScore3(ANSWER);
								break;
							default:
								temp.setScore3(Double.toString(MAX));
								break;
						}
						if ((int) temp.getScore2() == peak) {
							bedGraphMath.setSmooth(20, temp, overlaps);
							temp.setID("Peak_" + counter);
							if (i > 0) {
								temp.setUpstream(hmmr.get(i - 1));
							} else {
								temp.setUpstream(hmmr.get(i));
							}
							if (i < hmmr.size() - 1) {
								temp.setDownstream(hmmr.get(i + 1));
							} else {
								temp.setDownstream(hmmr.get(i));
							}
							counter++;
						}

					}
					/*
					 * report the bedgraph, is desired
					 */
					if (bg) {
						if (!BGScore) {
							bedgraph.println(temp.toString2());
						} else {
							bedgraph.println(temp.toString_ScoredBdg());
						}
					}
					/*
					 * report the peaks and summits, if desired
					 */
					if (peaks && (int) temp.getScore2() == peak
							&& temp.getLength() >= minLength && 
							Double.parseDouble(temp.getScore3()) >= threshold) {
						if (temp.getSummit() != null) {
							summits.println(temp.toString_ScoredSummit());
						}
						pks.println(temp.toString_gappedPeak());
					}

				}
			}
			
		}//for loop through chroms
		if (bg){
			bedgraph.close();
		}
		
		
		if (peaks){
			counter=0;

			for (TagNode tagNode : addBack) {
				String chrom = tagNode.getChrom();
				int start = tagNode.getStart();
				int stop = tagNode.getStop();

				pks.println(chrom + "\t" + start + "\t" + stop + "\t" + "HighCoveragePeak_" + counter + "\t.\t.\t0\t0\t255,0,0\t1\t" +
						tagNode.getLength() + "\t0\t-1\t-1\t-1");
			}
			
			pks.close();
			summits.close();
		}

		long endTime = System.currentTimeMillis();
		long total = (endTime - startTime) / 1000;
		log.println("Total time (seconds)= \t"+total);
		log.close();
	}//main
}

