package HMMR_ATAC;

import java.io.File;

public class ArgParser {
	
	private String[] args;
	
	//Required inputs
		private File bam = null;
		private File index = null;
		private  String genomeFile = null;
		private  String bigWig = null;
		
		//Optional Inputs
		private  String means;//comma separted list of intial mean values for frag dist
		private  String stddevs;//comma separted list of initial standard deviations for frag dist
		private  boolean fragEM = true; //whether or not to perform fragment dist em
		private  int minMapQ = 30; //minimum mapping quality of reads to keep
		private  int lower = 2; //lower bound for fold change range for choosing training sites
		private  int upper = 10; //upper bound for fold change range for choosing training sites
		private  int zscore = 100; //zscored read coverage to exclude from viterbi decoding
		private  String output; //output name 
		private  boolean peaks = false; // whether to print peaks
		private  boolean bg = true; // whether to print bedgraph
		private  String blacklist = null;
		private boolean keepDups = false;
		private int minLength = 200;
		private String scoreSys = "max";
		private String bgScore = "false";
		private int k = 4;
		
		public ArgParser(String[] a){
			args = a;
			set();
		}
		public String getBGScore(){return bgScore;}
		public String getScore(){return scoreSys;}
		public File getBam(){return bam;}
		public File getIndex(){return index;}
		public String getGenome(){return genomeFile;}
		public String getBigWig(){return bigWig;}
		public String getMeans(){return means;}
		public String getStd(){return stddevs;}
		public boolean getEM(){return fragEM;}
		public int getMinQ(){return minMapQ;}
		public int getLower(){return lower;}
		public int getUpper(){return upper;}
		public int getZscore(){return zscore;}
		public String getOutput(){return output;}
		public boolean getPeaks(){return peaks;}
		public boolean getBedgraph(){return bg;}
		public String getBlacklist(){return blacklist;}
		public boolean getKeepDups(){return keepDups;}
		public int getMinLength(){return minLength;}
		public int getK(){return k;}
		private void set(){
			for (int i = 0; i < args.length; i++) {

				switch (args[i].charAt((1))) {
				
				case'b':
					bam = new File(args[i+1]);
					i++;
					break;
				case'i':
					index = new File(args[i+1]);
					i++;
					break;
				case'g':
					genomeFile = args[i+1];
					i++;
					break;
				case'w':
					bigWig = args[i+1];
					i++;
					break;
				case'm':
					means = args[i+1];
					i++;
					break;
				case's':
					stddevs = args[i+1];
					i++;
					break;
				case'f':
					String temp = args[i+1].toLowerCase();
					if (temp.contains("t")){
						fragEM = true;
					}
					else{fragEM = false;}
					i++;
					break;
				case'q':
					minMapQ = Integer.parseInt(args[i+1]);
					i++;
					break;
				case'u':
					upper = Integer.parseInt(args[i+1]);
					i++;
					break;
				case'l':
					lower = Integer.parseInt(args[i+1]);
					i++;
					break;
				case'z':
					zscore = Integer.parseInt(args[i+1]);
					i++;
					break;
				case'o':
					output = args[i+1];
					i++;
					break;
				case'e':
					blacklist = args[i+1];
					i++;
					break;
				case'p':
					String t = args[i+1].toLowerCase();
					if (t.contains("t")){
						peaks=true;
					}
					else{peaks=false;}
					i++;
					break;
				case'k':
					k = Integer.parseInt(args[i+1]);
					i++;
					break;
				
				case'h':
					printUsage();
					//System.exit(0);
				case'-':
					switch(args[i].substring(2)){
					case "bam":
						bam = new File(args[i+1]);
						i++;
						break;
					case"index":
						index = new File(args[i+1]);
						i++;
						break;
					case"genome":
						genomeFile = args[i+1];
						i++;
						break;
					case"wig":
						bigWig = args[i+1];
						i++;
						break;
					case"means":
						means = args[i+1];
						i++;
						break;
					case"stddev":
						stddevs = args[i+1];
						i++;
						break;
					case"fragem":
						String t1 = args[i+1].toLowerCase();
						if(t1.contains("t")){
							fragEM=true;
						}
						else{fragEM = false;}
						i++;
						break;
					case"minmapq":
						minMapQ = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"upper":
						upper = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"lower":
						lower = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"zscore":
						zscore = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"output":
						output = args[i+1];
						i++;
						break;
					case"blacklist":
						blacklist = args[i+1];
						i++;
						break;
					case"peaks":
						String t2 = args[i+1].toLowerCase();
						if (t2.contains("t")){
							peaks=true;
						}
						else{peaks=false;}
						i++;
						break;
					case"bedgraph":
						String t3 = args[i+1].toLowerCase();
						if(t3.contains("t")){
							bg=true;
						}
						else{bg=false;}
						i++;
						break;
					case"minlen":
						minLength = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"score":
						scoreSys = args[i+1];
						i++;
						break;
					case"bgscore":
						bgScore = args[i+1];
						i++;
						break;
					case"kmeans":
						k = Integer.parseInt(args[i+1]);
						i++;
						break;
					case"help":
						printUsage();
						//System.exit(0);
					}
				}
			}
		}
		
		public void printUsage(){
			System.out.println("Usage: java -jar HMMR.jar");
			System.out.println("\nRequired Parameters:");
			System.out.println("\t-b , --bam <BAM> Sorted BAM file containing the ATAC-seq reads");
			System.out.println("\t-i , --index <BAI> Index file for the sorted BAM File");
			System.out.println("\t-g , --genome <GenomeFile> Two column, tab delimited file containing genome size stats");
			System.out.println("\t-w , --wig <BigWig> Whole genome big wig file created using all reads");
			System.out.println("\nOptional Parameters:");
			System.out.println("\t-m , --means <double> Comma separated list of initial mean values for the fragment distribution. Default = 50,200,400,600");
			System.out.println("\t-s , --stddev <double> Comma separated list of initial standard deviation values for fragment distribution. Default = 20,20,20,20");
			System.out.println("\t-f , --fragem <true/false> Whether to perform EM training on the fragment distribution. Default = True");
			System.out.println("\t-q , --minmapq <int> Minimum mapping quality of reads to keep. Default = 30");
			System.out.println("\t-u , --upper <int> Upper limit on fold change range for choosing training sites. Default = 10");
			System.out.println("\t-l , --lower <int> Lower limit on fold change range for choosing training sites. Default = 2");
			System.out.println("\t-z , --zscore <int> Zscored read depth to mask during Viterbi decoding. Default = 2 * upper FC value");
			System.out.println("\t-o , --output <File> Name for output files. Default = NA");
			System.out.println("\t-e , --blacklist <BED_File> bed file of blacklisted regions to exclude");
			System.out.println("\t-p , --peaks <true/false> Whether to report peaks in bed format. Default = false");
			System.out.println("\t-k , --kmeans <int> Number of States in the model. Default = 4. If not k=4, recommend NOT calling peaks, use bedgraph");
			System.out.println("\t--bedgraph <true/false> Whether to report whole genome bedgraph of all state anntations. Default = True");
			System.out.println("\t--minlen <int> Minimum length of open region to call peak. Note: -p , --peaks must be set. Default = 200");
			System.out.println("\t--score <max | ave | med | all> What type of score system to use for peaks. Can be used for ranking peaks. Default = max");
			System.out.println("\t--bgscore <true | false> Whether to add the HMMR score to each state annotation in bedgraph. Note: this adds considerable time. Default = False");
			System.out.println("\t-h , --help Print this help message and exit.");
			System.exit(0);
		}
		public static void main(String [] args){
			ArgParser a = new ArgParser(args);
			System.out.println(a.getPeaks());
			System.out.println(a.getEM());
			System.out.println(a.getBigWig());
		}
}
