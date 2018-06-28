package HMMR_ATAC;

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
import Node.TagNode;
import RobustHMM.KMeansToHMM;
import RobustHMM.RobustHMM;
import WigMath.FoldChange;
import WigMath.PullWigAboveCutoff;

public class Main_HMMR_Driver {

	//Required inputs
	private static File bam = null;
	private static File index = null;
	private static String genomeFile = null;
	private static String bigWig = null;
	
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
	private static int k = 4;
	private static int trim = 0;
	private static int vitWindow = 25000000;
	private static File modelFile;
	private static boolean stopAfterModel = false;
	
	private static String trainingRegions;
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		
		ArgParser p = new ArgParser(args);
		bam = p.getBam();
		index = p.getIndex();
		genomeFile=p.getGenome();
		bigWig=p.getBigWig();
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
		//For run time calculation
		Long startTime = System.currentTimeMillis();
		
		//Declare output name
		if (output==null){
			output = "NA";
		}
		PrintStream log = new PrintStream(output+".log");
		
		//Exit program if BAM file or Index file not given
		if (bam == null || index == null || genomeFile == null || bigWig == null){
			p.printUsage();
			System.exit(1);
		}
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
		
		//TODO: write method for generating whole genome read coverage
		
		/*
		 * Below is method to use BigWig input file
		 * Find regions with fold change within determined range to use as training sites.
		 * Find regions with zscore values above certain cutoff to exclude from viterbi.
		 */
		Hmm<ObservationVector> hmm=null;
		
		//
		FoldChange fc = new FoldChange(bigWig,upper,lower,genomeStats);
		PullWigAboveCutoff z = new PullWigAboveCutoff(bigWig,zscore,genomeStats);
		double genomeMean = fc.getMean();
		double genomeStd = fc.getSTD();
		
		ArrayList<TagNode> train = new MergeBed(fc.getResults()).getResults();
		
		ArrayList<TagNode> newTrain = new ArrayList<TagNode>();
		int maxTrain;
		if (train.size() > 1000){
			maxTrain = 1000;
		}
		else{
			maxTrain = train.size();
		}
		for (int i = 0;i < maxTrain;i++){
			newTrain.add(train.get(i));
		}
		train = newTrain;
		train = new ExtendBed(train,5000).getResults();
		fc = null;
		ArrayList<TagNode> exclude = new MergeBed(z.getResults()).getResults();
		ArrayList<TagNode> addBack = exclude;
		z = null;
		
		if (blacklist != null){
			exclude.addAll(black);
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
		
		FragPileupGen gen = new FragPileupGen(bam, index, train, mode, fragMeans, fragStddevs,minMapQ);
		TrackHolder holder = new TrackHolder(gen.transformTracks(gen.getAverageTracks()),trim);// 8/30/16 transformed tracks
		
		gen = null;
		
		log.println("Training Fragment Pileup completed");
		
		/*
		 * Create the initial model using KMeans and then refine it using Baum-Welch
		 */
		
		KMeansToHMM kmeans = new KMeansToHMM(holder.getDataSet(),k,Integer.MAX_VALUE,true,true,true);
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
		ArrayList<TagNode> genomeAnnotation = new ArrayList<TagNode>();
		for (int i = 0;i < vitBed.size();i++){
			if (vitBed.get(i).getLength() >= 10){
				ArrayList<TagNode> tempBed = new ArrayList<TagNode>();
				tempBed.add(vitBed.get(i));
				FragPileupGen vGen = new FragPileupGen(bam, index, tempBed, mode, fragMeans, fragStddevs,minMapQ);
				TrackHolder vHolder = new TrackHolder(vGen.transformTracks(vGen.getAverageTracks()),trim);// 8/30/16 transformed tracks
				//Reverse the tracks as a test 11/15/17
				vGen = null;
			
				RobustHMM HMM = new RobustHMM(vHolder.getObs(),null,hmm,false,0,"Vector",0);
				int[] states = HMM.getStates();
				//Reverse the states for creating peaks
				//ArrayUtils.reverse(states);
				
				int start = vitBed.get(i).getStart();
				
				ArrayList<PileupNode2> pile = new ArrayList<PileupNode2>();
				for (int a = 0;a < states.length;a++){
					PileupNode2 pNode = new PileupNode2(start+(a*10),(double)states[a],vitBed.get(i).getChrom());
					pile.add(pNode);
					
					
					//out.println(pNode.getChrom()+"\t"+pNode.getBase()+"\t"+(pNode.getBase()+10)+"\t"+"E"+(int)pNode.getScore());
				}
				genomeAnnotation.addAll(new PileupToBedGraph(pile,10).getBedGraph());
				
				log.println(i+" round viterbi done");
			}
		}
		//out.close();
		
		/*
		 * Print Genome-wide bedgraph of state annotations if desired
		 */
		
		if (bg){
			PrintStream bedgraph = new PrintStream(output+".bedgraph");
		
			for (int i = 0;i < genomeAnnotation.size();i++){
				String chr = genomeAnnotation.get(i).getChrom();
				int start = genomeAnnotation.get(i).getStart();
				int stop = genomeAnnotation.get(i).getStop();
				double score = genomeAnnotation.get(i).getScore2();
				if (BGScore){
					GetSignal sig = new GetSignal(bigWig,chr,start,stop);
					bedgraph.println(chr+"\t"+start+"\t"+stop+"\t"+"E"+(int)score+"\t"+sig.getMax());
				}
				else{
					bedgraph.println(chr+"\t"+start+"\t"+stop+"\t"+"E"+(int)score);
				}
			}
			bedgraph.close();
		}
		
		/*
		 * Print Peak file and summit file if desired
		 */
		
		if (peaks){
			PrintStream pks = new PrintStream(output+"_peaks.gappedPeak");
			@SuppressWarnings("resource")
			PrintStream summits = new PrintStream(output+"_summits.bed");
			pks.println("track name="+output+"_peaks.gappedPeak type=gappedPeak");
			int counter=0;
			for (int i = 1; i < genomeAnnotation.size()-1;i++){
				if (genomeAnnotation.get(i).getLength() >= minLength){
					String chr = genomeAnnotation.get(i).getChrom();
					int start = genomeAnnotation.get(i).getStart();
					int stop = genomeAnnotation.get(i).getStop();
					int score = (int)genomeAnnotation.get(i).getScore2();
					if (chr.equals(genomeAnnotation.get(i-1).getChrom()) && chr.equals(genomeAnnotation.get(i+1).getChrom())){
						if (score == peak){
							if(genomeAnnotation.get(i-1).getStop() == start && 
									genomeAnnotation.get(i+1).getStart() == stop){
								counter+=1;
								GetSignal sig = new GetSignal(bigWig,chr,start-60,stop+60);
								String value = Double.toString(sig.getMax());
								if (scoreSys.equals("ave")){
									value = Double.toString(sig.getScore());
								}
								else if (scoreSys.equals("fc")){
									value = Double.toString(sig.getScore()/genomeMean);
								}
								else if (scoreSys.equals("zscore")){
									value = Double.toString((sig.getScore()-genomeMean)/genomeStd);
								}
								else if (scoreSys.equals("med")){
									value = Double.toString(sig.getMedian());
								}
								else if (scoreSys.equals("all")){
									value = Double.toString(sig.getScore())+"_"+
											Double.toString(sig.getMax())+"_"+Double.toString(sig.getMedian())+
											"_"+Double.toString(sig.getScore()/genomeMean)+"_"+
											Double.toString((sig.getScore()-genomeMean)/genomeStd);
								}
								
								// Find summits
								//New Way 1/24/17 - Report Position where gaussian-smoothed max score is found
								if (sig.getMaxPos() != null){
									
									sig.findSmooth(20);
									TagNode best = sig.getMaxPos();
									summits.println(best.toString()+"\t"+"Peak_"+counter+"\t"+value);
								}
								//End report summits
								
								int peakStart=genomeAnnotation.get(i-1).getStart();
								int peakStop= genomeAnnotation.get(i+1).getStop();
								String middleValues = "3"+"\t"+
									"1"+","+genomeAnnotation.get(i).getLength()+","+
										"1"+"\t"+"0,"+
									(genomeAnnotation.get(i).getStart()-genomeAnnotation.get(i-1).getStart())+","+
										((genomeAnnotation.get(i+1).getStop()-genomeAnnotation.get(i-1).getStart())-1);
								/*
								if (genomeAnnotation.get(i-1).getScore2()== 0 && genomeAnnotation.get(i+1).getScore2()!=0){
									//non-standard peak, only downstream nuc
									peakStart = genomeAnnotation.get(i).getStart();
									peakStop = genomeAnnotation.get(i+1).getStop();
									middleValues = "2\t"+genomeAnnotation.get(i).getLength()+","+"1"+"\t"+"0"+","+((genomeAnnotation.get(i+1).getStop()-genomeAnnotation.get(i).getStart())-1);
								}
								else if (genomeAnnotation.get(i-1).getScore2()!= 0 && genomeAnnotation.get(i+1).getScore2()==0){
									//non-standard peak, only upstream nuc
									peakStart = genomeAnnotation.get(i-1).getStart();
									peakStop = genomeAnnotation.get(i).getStop();
									middleValues = "2\t"+"1"+","+genomeAnnotation.get(i).getLength()+"\t"+"0"+","+(genomeAnnotation.get(i).getStart()-genomeAnnotation.get(i-1).getStart());

								}
								else if (genomeAnnotation.get(i-1).getScore2()== 0 && genomeAnnotation.get(i+1).getScore2()==0){
									//non-standard peak, no nucs
									peakStart = genomeAnnotation.get(i).getStart();
									peakStop = genomeAnnotation.get(i).getStop();
									middleValues = "1\t"+genomeAnnotation.get(i).getLength()+"\t"+"0";
								}
								*/
								pks.println(chr+"\t"+peakStart+"\t"+
									peakStop+"\t"+"Peak_"+counter+"\t"+".\t."+"\t"+genomeAnnotation.get(i).getStart()+"\t"+
										genomeAnnotation.get(i).getStop()+"\t"+"255,0,0"+"\t"+middleValues+"\t"+value+"\t"+"-1\t-1");
							}
						}
					}
				}
			}
			//add back the high coverage regions
			counter=0;
			for (int i = 0;i < addBack.size();i++){
				String chrom = addBack.get(i).getChrom();
				int start = addBack.get(i).getStart();
				int stop = addBack.get(i).getStop();
				GetSignal sig = new GetSignal(bigWig,chrom,start,stop);
				String value = Double.toString(sig.getMax());
				counter+=1;
				if (scoreSys.equals("ave")){
					value = Double.toString(sig.getScore());
				}
				else if (scoreSys.equals("med")){
					value = Double.toString(sig.getMedian());
				}
				else if (scoreSys.equals("all")){
					value = sig.getScore()+"_"+sig.getMax()+"_"+sig.getMedian();
				}
				pks.println(chrom+"\t"+start+"\t"+stop+"\t"+"HighCoveragePeak_"+counter+"\t.\t.\t0\t0\t255,0,0\t1\t"+
				addBack.get(i).getLength()+"\t0\t"+value+"\t-1\t-1");
			}
			pks.close();
		}
		Long endTime = System.currentTimeMillis();
		Long total = (endTime - startTime) / 1000;
		log.println("Total time (seconds)= \t"+total);
		log.close();
	}//main
	
	
}

