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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

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

	//Required inputs
	private static File bam = null;
	private static File index = null;
	private static String genomeFile = null;
	//private static String bigWig = null;//10_13_18
	
	//Optional Inputs
	private static String means;//comma separted list of intial mean values for frag dist
	private static String stddevs;//comma separted list of initial standard deviations for frag dist
	private static boolean fragEM = true; //whether or not to perform fragment dist em
	private static int minMapQ = 30; //minimum mapping quality of reads to keep
	
	private static int lower = 10; //lower bound for fold change range for choosing training sites
	private static int upper = 20; //upper bound for fold change range for choosing training sites
	private static int zscore = 100; //zscored read coverage to exclude from viterbi decoding
	private static String output; //output name 
	private static boolean peaks; // whether to print peaks
	private static boolean bg; // whether to print bedgraph
	private static String blacklist = null;
	private static int minLength;
	private static String scoreSys;
	private static boolean BGScore;
	private static int k = 3;
	private static int trim = 0;
	private static int vitWindow = 25000000;
	private static File modelFile;
	private static boolean stopAfterModel = false;
	private static boolean printHMMRTracks = false;
	private static int maxTrain = 1000;
	private static boolean rmDup = true;
	private static boolean printExclude = false;
	private static boolean printTrain = true;
	private static long randomTrainSeed=10151;
	private static double threshold=0;
	
	private static String trainingRegions;
	
	/*
	 * Version number. Change as needed
	 */
	private static String versionNum = "1.2.10";
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		
		ArgParser p = new ArgParser(args,versionNum);
		bam = p.getBam();
		index = p.getIndex();
		genomeFile=p.getGenome();
		//bigWig=p.getBigWig();//10_13_18
		means = p.getMeans();
		stddevs=p.getStd();
		fragEM=p.getEM();
		minMapQ=p.getMinQ();
		lower=p.getLower();
		upper=p.getUpper();
		zscore=p.getZscore();
		output=p.getOutput();
		peaks=p.getPeaks();
		bg=p.getBedgraph();
		blacklist=p.getBlacklist();
		minLength = p.getMinLength();
		scoreSys = p.getScore();
		BGScore = p.getBGScore();
		k = p.getK();
		trim=p.getTrim();
		trainingRegions = p.getTrainingRegions();
		vitWindow = p.getWindow();
		modelFile = p.getModelFile();
		maxTrain = p.getMaxTrain();
		rmDup = p.getRemoveDuplicates();
		printExclude = p.getPrintExclude();
		printTrain = p.getPrintTrain();
		randomTrainSeed = p.getRandomTrainSeed();
		threshold=p.getThreshold();
//		printHMMRTracks = p.getPrintHMMRTracks(); 
		//For run time calculation
		Long startTime = System.currentTimeMillis();
		
		//Declare output name
		if (output==null){
			output = "NA";
		}
		PrintStream log = new PrintStream(output+".log");
		
		//Exit program if BAM file or Index file not given
		if (bam == null || index == null || genomeFile == null ){//|| bigWig == null){//10_13_18
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
		gReader = null;
		
		//Read in blacklisted if inputted
		ArrayList<TagNode> black = null;
		if (blacklist != null){
			black = new bedFileReader(blacklist).getData();
			
		}
		
		/*
		 * Set fragment length distribution parameters. 
		 * Use inputed values to set initial values, if provided. 
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
			for (int i = 0;i < fragStddevs.length;i++){
				fragStddevs[i] = 20.0;
			}
		}
		
		
		/*
		 * Pull the lengths from the read data. Fragment distribution uses fragments with length > 100 to train the 
		 nucleosome distributions. the short distribution is set to 50 and remains unchanged
		 Only occurs if EM training occurs, else use the default settings
		*/
		
		
		if (fragEM){
			pullLargeLengths puller = new pullLargeLengths(bam, index, minMapQ, genomeStats,java.util.Arrays.copyOf(fragMeans, 4));
			double[] lengths = puller.getSampledLengths(10);
			double[] weights = puller.getWeights();
			puller = null;
			
			//Perform EM training
			
			HMMR_EM em = new HMMR_EM(weights,java.util.Arrays.copyOfRange(fragMeans, 1,4),
					java.util.Arrays.copyOfRange(fragStddevs, 1, 4),lengths);
			em.learn();
			double[] tempMeans = em.getMeans();
			double[] tempLam = em.getLamda();
			em = null;
			for (int i = 0;i < tempMeans.length;i++){
				//This will update the parameters IFF they were updated. If they become NaN, leave as default
				if(!Double.isNaN(tempMeans[i]) && !Double.isNaN(tempLam[i])){
					fragMeans[i+1] = tempMeans[i];
					fragStddevs[i+1] = tempLam[i];
				}
			}
			em = null;lengths = null;tempMeans=null;tempLam=null;
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
		 * NOTE: this method does not currently exist. Instead we use an inputed big wig file to find fold change 
		 * training regions and zscored masked regions
		 */
		
		
		
		/*
		 * Below is method to use BigWig input file
		 * Find regions with fold change within determined range to use as training sites.
		 * Find regions with zscore values above certain cutoff to exclude from viterbi.
		 */
		Hmm<ObservationVector> hmm=null;
		
		// Comment out 10_13_18 to eliminate bigwig requirement
		//FoldChange fc = new FoldChange(bigWig,upper,lower,genomeStats);
		//PullWigAboveCutoff z = new PullWigAboveCutoff(bigWig,zscore,genomeStats);
		
		pileup pileupData = new pileup(new SplitBed(genomeStats,vitWindow).getResult(), 0, bam, index, 0,rmDup);
//		pileup pileupData = new pileup(genomeStats, 0, bam, index, minMapQ);
		bedGraphMath fc = new bedGraphMath(pileupData.getBedGraph());
		
		//calculate the cpm scaling factor for input into FragPileupGen. Use cpmScale=1 for no scaling
		double cpmScale = pileupData.getCPMScale()/1000000;
		log.println("ScalingFactor\t"+cpmScale);
		
		pileupData = null;
		
		
		//ArrayList<TagNode> bdg2 = fc.getBedGraph();
		//for (int i = 0; i < bdg2.size();i++){
		//	System.out.println(bdg2.get(i).toString3()+"\tBedgraph Signal");
		//}
		
		
		double genomeMean = fc.getMean();
		double genomeStd = fc.getSTD();
		
		ArrayList<TagNode> train = new MergeBed(fc.getBetweenRanges(upper, lower)).getResults();
		
		
		
		ArrayList<TagNode> newTrain = new ArrayList<TagNode>();
//		int maxTrain;
		if (train.size() < maxTrain){
			maxTrain = train.size();
		}
//		else{
//			maxTrain = train.size();
//		}
		
		//Shuffle training list before choosing.
		Collections.shuffle(train, new Random(randomTrainSeed));
		for (int i = 0;i < maxTrain;i++){
			newTrain.add(train.get(i));
		}
		
		
		
		train = newTrain;
		train = new ExtendBed(train,5000).getResults();
		//fc = null; 10_13_18
		ArrayList<TagNode> exclude = new MergeBed(fc.getAboveZscore(zscore)).getResults();
		ArrayList<TagNode> addBack = exclude;
		//z=null; 10_13_18
		
		if (blacklist != null){
			exclude.addAll(black);
			exclude = new MergeBed(exclude).getResults();//added 8-3-18
		}
			
		if(printExclude){
			PrintStream ex = new PrintStream(output+"_excluded.bed");
			for (int c = 0;c < exclude.size();c++){
				ex.println(exclude.get(c).toString()+"\t"+"exclude");
			}
			ex.close();
		}
		
		
		newTrain = new ArrayList<TagNode>();
		for (int i = 0;i < train.size();i++){
			int counter = 0;
			TagNode node1 = train.get(i);
			for (int a = 0;a < exclude.size();a++){
				TagNode node2 = exclude.get(a);
				if (SubtractBed.overlap(node1, node2).hasHit()){
					counter++;
				}
			}
			if (counter == 0){
				newTrain.add(node1);
					
			}
		}
		
		train = newTrain;
		
		//Below added on 12/13/16
		// Allows user to use training set for model generation
		if (trainingRegions != null){
			bedFileReader trainReader = new bedFileReader(trainingRegions);
			train = trainReader.getData();
		}
		//Above added 12/13/16
		
		log.println("Training Regions found and Zscore regions for exclusion found");
		
		
		/*
		 * Create the fragment pileup tracks using the training set and the fragment distribution parameters
		 */
		if (modelFile == null){
		
			
			if (printTrain) {
				PrintStream tr = new PrintStream(output + "_training.bed");
				for (int i = 0; i < train.size(); i++) {
					tr.println(train.get(i).toString() + "\t" + "training");
				}
				tr.close();
			}
			FragPileupGen gen = new FragPileupGen(bam, index, train, mode, fragMeans, fragStddevs,minMapQ,rmDup,cpmScale);
			TrackHolder holder = new TrackHolder((gen.transformTracks(gen.scaleTracks(gen.getAverageTracks()))),trim);// 8/30/16 transformed tracks 
			//	7/16/18 transformation removed after testing showed it has little effect with new weighted procedure 
			
			
			gen = null;
			
			log.println("Training Fragment Pileup completed");
			
			/*
			 * Create the initial model using KMeans and then refine it using Baum-Welch
			 */
			
			KMeansToHMM kmeans = new KMeansToHMM(holder.getDataSet(),k,Integer.MAX_VALUE,true,true,true);
			log.println("Kmeans Model:\n"+kmeans.getHMM().toString()); // added 7-13-18
			//System.out.println(kmeans.getHMM().toString());
			
			hmm = new BaumWelch(kmeans.getHMM(),holder.getBWObs(),1000).build();
		
			//System.out.println(hmm.toString());
			kmeans = null; holder=null;
		
		}
		/*
		 * Use input model if available
		 */
		
		if (modelFile != null){
			hmm = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(modelFile));
		}

		
		/*
		 * Identify peak state as the state with the highest short read signal.
		 * Identify flanking nucleosome state as the state with the second highest mono read signal.
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
		File outputModel = new File(output+".model");
		FileOutputStream outModel = new FileOutputStream(outputModel);
		HmmBinaryWriter.write(outModel, hmm);
		outModel.close();
		log.println("Model created and refined. See "+output+".model");
		log.println("Model:\n"+hmm.toString());
		
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
		
		ArrayList<TagNode> split = new SplitBed(genomeStats,vitWindow).getResult();
		genomeStats = null;
		
		
		
		/*
		 * Subtract excluded regions from the split genome for Viterbi
		 */
		
		ArrayList<TagNode> vitBed = new SubtractBed(split,exclude).getResults();
		split=null;exclude=null;
		
		log.println("Genome split and subtracted masked regions");
		
		/*
		 * Run viterbi on the whole genome
		 */
		//PrintStream out = new PrintStream(output+".pileup");
		PrintStream NFR = null;
		PrintStream MONO = null;
		PrintStream DI=null;
		PrintStream TRI=null;
		if(printHMMRTracks){
			 NFR = new PrintStream(output+"_nfr.bedgraph");
			 MONO = new PrintStream(output+"_mono.bedgraph");
			 DI = new PrintStream(output+"_di.bedgraph");
			 TRI = new PrintStream(output+"_tri.bedgraph");
		}
		ArrayList<TagNode> genomeAnnotation = new ArrayList<TagNode>();
		for (int i = 0;i < vitBed.size();i++){
			if (vitBed.get(i).getLength() >= 10){
				ArrayList<TagNode> tempBed = new ArrayList<TagNode>();
				tempBed.add(vitBed.get(i));
				FragPileupGen vGen = new FragPileupGen(bam, index, tempBed, mode, fragMeans, fragStddevs,minMapQ,rmDup,cpmScale);
				TrackHolder vHolder = new TrackHolder(vGen.transformTracks(vGen.scaleTracks(vGen.getAverageTracks())),trim);// 8/30/16 transformed tracks
				
				if (printHMMRTracks){
					HMMRTracksToBedgraph tracks = new HMMRTracksToBedgraph(vHolder.getRawData(),vitBed.get(i),10);
					ArrayList<TagNode> nfr = tracks.getShort();
					ArrayList<TagNode> mono = tracks.getMono();
					ArrayList<TagNode> di = tracks.getDi();
					ArrayList<TagNode> tri = tracks.getTri();
					
					
					if (nfr != null) {
						for (int w = 0; w < nfr.size(); w++) {
							NFR.println(nfr.get(w).toString2());
						}
					}
					if (mono != null) {
						for (int d = 0; d < mono.size(); d++) {
							MONO.println(mono.get(d).toString2());
						}
					}
					if (di != null) {
						for (int e = 0; e < di.size(); e++) {
							DI.println(di.get(e).toString2());
						}
					}
					if (tri != null) {
						for (int f = 0; f < tri.size(); f++) {
							TRI.println(tri.get(f).toString2());
						}
					}
					
				}
				vGen = null;
			
				RobustHMM HMM = new RobustHMM(vHolder.getObs(),null,hmm,false,0,"Vector",0);
				int[] states = HMM.getStates();
				//Reverse the states for creating peaks
				//ArrayUtils.reverse(states);
				
				int start = vitBed.get(i).getStart();
				int remainder = vitBed.get(i).getLength() % 10;
				ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
				int a;
				for (a = 0;a < states.length-1;a++){
					//PileupNode2 pNode = new PileupNode2(start+(a*10),(double)states[a],vitBed.get(i).getChrom());
					PileupNode2 pNode = new PileupNode2(start+(a*10),(double)states[a],vitBed.get(i).getChrom());
					pile.add(pNode);
					
					
					//out.println(pNode.getChrom()+"\t"+pNode.getBase()+"\t"+(pNode.getBase()+10)+"\t"+"E"+(int)pNode.getScore());
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
		/**
		 * Report the final results as peaks, bedgraphs and summits, if desired
		 */
		PrintStream bedgraph=null;
		if (bg){
			 bedgraph = new PrintStream(output+".bedgraph");
		}
		PrintStream pks=null;
		PrintStream summits=null;
		if (peaks){
			 pks = new PrintStream(output+"_peaks.gappedPeak");
			 summits = new PrintStream(output+"_summits.bed");
		}
		HashMap<String,ArrayList<TagNode>> bdg = fc.getMappedBedgraph();
//		HashMap<String,ArrayList<TagNode>> bdg = null;
		
		fc=null;
		HashMap<String,ArrayList<TagNode>> hmmrBdg = bedGraphMath.toMap(genomeAnnotation);
		int counter=1;
		for (String chr : hmmrBdg.keySet()){
			ArrayList<TagNode> hmmr = hmmrBdg.get(chr);
			ArrayList<TagNode> signal = bdg.get(chr);
			if (signal != null) {
				Collections.sort(hmmr, TagNode.basepairComparator);
				
				Collections.sort(signal, TagNode.basepairComparator);
				
				int index = 0;
				for (int i = 0; i < hmmr.size(); i++) {
					TagNode temp = hmmr.get(i);

					/**
					 * Execute the scoring commands if the state is a peak or if bedgraph scoring is on
					 */
					if ((int) temp.getScore2() == peak || BGScore) {
						boolean hasHadOverlap = false;
						ArrayList<TagNode> overlaps = new ArrayList<TagNode>();
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
						double ZSCORE = (scores.getMean() - genomeMean)
								/ genomeStd;
						double FOLDCHANGE = scores.getMean() / genomeMean;
						
						if (scoreSys.equals("ave")) {
							temp.setScore3(Double.toString(MEAN));
						} else if (scoreSys.equals("fc")) {
							temp.setScore3(Double.toString(FOLDCHANGE));
						} else if (scoreSys.equals("zscore")) {
							temp.setScore3(Double.toString(ZSCORE));
						} else if (scoreSys.equals("med")) {
							temp.setScore3(Double.toString(MEDIAN));
						} else if (scoreSys.equals("all")){
							String ANSWER = MAX+"_"+MEAN+"_"+MEDIAN+"_"+ZSCORE+"_"+FOLDCHANGE;
							temp.setScore3(ANSWER);
						} else {
							temp.setScore3(Double.toString(MAX));
						}
						if ((int) temp.getScore2() == peak) {
							temp = bedGraphMath.setSmooth(20, temp, overlaps);
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
					/**
					 * report the bedgraph, is desired
					 */
					if (bg) {
						if (!BGScore) {
							bedgraph.println(temp.toString2());
						} else {
							bedgraph.println(temp.toString_ScoredBdg());
						}
					}
					/**
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
			
			for (int i = 0;i < addBack.size();i++){
				String chrom = addBack.get(i).getChrom();
				int start = addBack.get(i).getStart();
				int stop = addBack.get(i).getStop();
				
				pks.println(chrom+"\t"+start+"\t"+stop+"\t"+"HighCoveragePeak_"+counter+"\t.\t.\t0\t0\t255,0,0\t1\t"+
				addBack.get(i).getLength()+"\t0\t-1\t-1\t-1");
			}
			
			pks.close();
			summits.close();
		}
		//old way
		
		/*
		 * Print Genome-wide bedgraph of state annotations if desired
		 */
//		if (bg){
//			PrintStream bedgraph = new PrintStream(output+".bedgraph");
//		
//			for (int i = 0;i < genomeAnnotation.size();i++){
//				String chr = genomeAnnotation.get(i).getChrom();
//				int start = genomeAnnotation.get(i).getStart();
//				int stop = genomeAnnotation.get(i).getStop();
//				double score = genomeAnnotation.get(i).getScore2();
//				if (BGScore){
//					//GetSignal sig = new GetSignal(bigWig,chr,start,stop);10_13_18
//					fc.set(genomeAnnotation.get(i));
//					bedgraph.println(chr+"\t"+start+"\t"+stop+"\t"+"E"+(int)score+"\t"+fc.getMax());//10_13_18
//					//bedgraph.println(chr+"\t"+start+"\t"+stop+"\t"+"E"+(int)score+"\t"+sig.getMax());10_13_18
//				}
//				else{
//					bedgraph.println(chr+"\t"+start+"\t"+stop+"\t"+"E"+(int)score);
//				}
//			}
//			bedgraph.close();
//		}
//		
//		/*
//		 * Print Peak file and summit file if desired
//		 */
//		
//		if (peaks){
//			PrintStream pks = new PrintStream(output+"_peaks.gappedPeak");
//			@SuppressWarnings("resource")
//			PrintStream summits = new PrintStream(output+"_summits.bed");
//			pks.println("track name="+output+"_peaks.gappedPeak type=gappedPeak");
//			int counter=0;
//			for (int i = 1; i < genomeAnnotation.size()-1;i++){
//				if (genomeAnnotation.get(i).getLength() >= minLength){
//					String chr = genomeAnnotation.get(i).getChrom();
//					int start = genomeAnnotation.get(i).getStart();
//					int stop = genomeAnnotation.get(i).getStop();
//					int score = (int)genomeAnnotation.get(i).getScore2();
//					if (chr.equals(genomeAnnotation.get(i-1).getChrom()) && chr.equals(genomeAnnotation.get(i+1).getChrom())){
//						if (score == peak){
//							if(genomeAnnotation.get(i-1).getStop() == start && 
//									genomeAnnotation.get(i+1).getStart() == stop){
//								counter+=1;
//								//GetSignal sig = new GetSignal(bigWig,chr,start-60,stop+60);10_13_18
//								fc.set(new TagNode(chr,start-60,stop+60));//10_13_18
//								//String value = Double.toString(sig.getMax());//10_13_18
//								String value = Double.toString(fc.getMax());//10_13_18
//								if (scoreSys.equals("ave")){
//									//value = Double.toString(sig.getScore());//10_13_18
//									value = Double.toString(fc.getScore());//10_13_18
//								}
//								else if (scoreSys.equals("fc")){
//									//value = Double.toString(sig.getScore()/genomeMean);//10_13_18
//									value = Double.toString(fc.getScore()/genomeMean);//10_13_18
//								}
//								else if (scoreSys.equals("zscore")){
//									//value = Double.toString((sig.getScore()-genomeMean)/genomeStd);//10_13_18
//									value = Double.toString((fc.getScore()-genomeMean)/genomeStd);//10_13_18
//								}
//								else if (scoreSys.equals("med")){
//									//value = Double.toString(sig.getMedian());//10_13_18
//									value = Double.toString(fc.getMedian());//10_13_18
//								}
//								else if (scoreSys.equals("all")){
//									//value = Double.toString(sig.getScore())+"_"+
//										//	Double.toString(sig.getMax())+"_"+Double.toString(sig.getMedian())+
//										//	"_"+Double.toString(sig.getScore()/genomeMean)+"_"+
//										//	Double.toString((sig.getScore()-genomeMean)/genomeStd);//10_13_18
//									value = Double.toString(fc.getScore())+"_"+
//											Double.toString(fc.getMax())+"_"+Double.toString(fc.getMedian())+
//											"_"+Double.toString(fc.getScore()/genomeMean)+"_"+
//											Double.toString((fc.getScore()-genomeMean)/genomeStd);//10_13_18
//								}
//								
//								// Find summits
//								//New Way 1/24/17 - Report Position where gaussian-smoothed max score is found
//								//if (sig.getMaxPos() != null){//10_13_18
//								if (fc.getMaxPos() != null){//10_13_18	
//									//sig.findSmooth(20);//10_13_18
//									fc.setSmooth(20, new TagNode(chr,start-60,stop+60));//10_13_18
//									
//									//TagNode best = sig.getMaxPos();//10_13_18
//									TagNode best = fc.getMaxPos();//10_13_18
//									summits.println(best.toString()+"\t"+"Peak_"+counter+"\t"+value);
//								}
//								//End report summits
//								
//								int peakStart=genomeAnnotation.get(i-1).getStart();
//								int peakStop= genomeAnnotation.get(i+1).getStop();
//								String middleValues = "3"+"\t"+
//									"1"+","+genomeAnnotation.get(i).getLength()+","+
//										"1"+"\t"+"0,"+
//									(genomeAnnotation.get(i).getStart()-genomeAnnotation.get(i-1).getStart())+","+
//										((genomeAnnotation.get(i+1).getStop()-genomeAnnotation.get(i-1).getStart())-1);
//								/*
//								if (genomeAnnotation.get(i-1).getScore2()== 0 && genomeAnnotation.get(i+1).getScore2()!=0){
//									//non-standard peak, only downstream nuc
//									peakStart = genomeAnnotation.get(i).getStart();
//									peakStop = genomeAnnotation.get(i+1).getStop();
//									middleValues = "2\t"+genomeAnnotation.get(i).getLength()+","+"1"+"\t"+"0"+","+((genomeAnnotation.get(i+1).getStop()-genomeAnnotation.get(i).getStart())-1);
//								}
//								else if (genomeAnnotation.get(i-1).getScore2()!= 0 && genomeAnnotation.get(i+1).getScore2()==0){
//									//non-standard peak, only upstream nuc
//									peakStart = genomeAnnotation.get(i-1).getStart();
//									peakStop = genomeAnnotation.get(i).getStop();
//									middleValues = "2\t"+"1"+","+genomeAnnotation.get(i).getLength()+"\t"+"0"+","+(genomeAnnotation.get(i).getStart()-genomeAnnotation.get(i-1).getStart());
//
//								}
//								else if (genomeAnnotation.get(i-1).getScore2()== 0 && genomeAnnotation.get(i+1).getScore2()==0){
//									//non-standard peak, no nucs
//									peakStart = genomeAnnotation.get(i).getStart();
//									peakStop = genomeAnnotation.get(i).getStop();
//									middleValues = "1\t"+genomeAnnotation.get(i).getLength()+"\t"+"0";
//								}
//								*/
//								pks.println(chr+"\t"+peakStart+"\t"+
//									peakStop+"\t"+"Peak_"+counter+"\t"+".\t."+"\t"+genomeAnnotation.get(i).getStart()+"\t"+
//										genomeAnnotation.get(i).getStop()+"\t"+"255,0,0"+"\t"+middleValues+"\t"+value+"\t"+"-1\t-1");
//							}
//						}
//					}
//				}
//			}
//			//add back the high coverage regions
//			counter=0;
//			for (int i = 0;i < addBack.size();i++){
//				String chrom = addBack.get(i).getChrom();
//				int start = addBack.get(i).getStart();
//				int stop = addBack.get(i).getStop();
//				//GetSignal sig = new GetSignal(bigWig,chrom,start,stop);//10_13_18
//				fc.set(new TagNode(chrom,start,stop));//10_13_18
//				//String value = Double.toString(sig.getMax());//10_13_18
//				String value = Double.toString(fc.getMax());//10_13_18
//				counter+=1;
//				if (scoreSys.equals("ave")){
//					//value = Double.toString(sig.getScore());//10_13_18
//					value = Double.toString(fc.getScore());//10_13_18
//				}
//				else if (scoreSys.equals("med")){
//					//value = Double.toString(sig.getMedian());//10_13_18
//					value = Double.toString(fc.getMedian());//10_13_18
//				}
//				else if (scoreSys.equals("all")){
//					//value = sig.getScore()+"_"+sig.getMax()+"_"+sig.getMedian();//10_13_18
//					value = fc.getScore()+"_"+fc.getMax()+"_"+fc.getMedian();//10_13_18
//				}
//				pks.println(chrom+"\t"+start+"\t"+stop+"\t"+"HighCoveragePeak_"+counter+"\t.\t.\t0\t0\t255,0,0\t1\t"+
//				addBack.get(i).getLength()+"\t0\t"+value+"\t-1\t-1");
//			}
//			pks.close();
//		}
//		
		Long endTime = System.currentTimeMillis();
		Long total = (endTime - startTime) / 1000;
		log.println("Total time (seconds)= \t"+total);
		log.close();
	}//main
	
	
}

