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

import java.io.File;
import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;


public class ArgParser {

	private static class FileConverter implements IStringConverter<File> {
		@Override
		public File convert(String value) {
			return new File(value);
		}
	}

	//Required inputs

	@Parameter(names={"-b", "--bam"}, converter=FileConverter.class, required = true, description="Sorted BAM file containing the ATAC-seq reads")
	public File bam = null;

	@Parameter(names={"-i", "--index"}, converter=FileConverter.class, required = true, description="Index file for the sorted BAM File")
	public File index = null;

	@Parameter(names={"-g", "--genome"}, converter=FileConverter.class, required = true, description="Two column, tab delimited file containing genome size stats")
	public File genomeFile = null;

	//Optional Inputs

	@Parameter(names={"-m", "--means"}, arity=4, description="Comma separated list of initial mean values for the fragment distribution. Default = 50,200,400,600")
	public double[] means = {50.0, 200.0, 400.0, 600.0};

	@Parameter(names={"-s", "--stddev"}, arity=4, description="Comma separated list of initial standard deviation values for fragment distribution. Default = 20,20,20,20")
	public double[] stddevs = {20.0, 20.0, 20.0, 20.0};

	@Parameter(names={"-f", "--fragem"}, arity=1, description="Whether to perform EM training on the fragment distribution. Default = True")
	public boolean fragEM = true;

	@Parameter(names={"-q", "--minmapq"}, description="Minimum mapping quality of reads to keep. Default = 30")
	public int minMapQ = 30;

	@Parameter(names={"-u", "--upper"}, description="Upper limit on fold change range for choosing training sites. Default = 20")
	public int lower = 10;

	@Parameter(names={"-l", "--lower"}, description="Lower limit on fold change range for choosing training sites. Default = 10")
	public int upper = 20;

	@Parameter(names={"-z", "--zscore"}, description="Zscored read depth to mask during Viterbi decoding. Default = 100")
	public int zscore = 100;

	@Parameter(names={"o", "--output"}, description="Name for output files. Default = NA")
	public String output = "NA";

	@Parameter(names={"-p", "--peaks"}, arity=1 , description="Whether to report peaks in bed format. Default = true")
	public boolean peaks = true;

	@Parameter(names={"--bedgraph"}, arity=1 , description="Whether to report whole genome bedgraph of all state annotations. Default = false")
	public boolean bg = false;

	@Parameter(names={"e", "--blacklist"}, converter=FileConverter.class, description="bed file of blacklisted regions to exclude")
	public File blacklist = null;

	@Parameter(names={"--minlen"}, description="Minimum length of open region to call peak. Note: -p , --peaks must be set. Default = 200")
	public int minLength = 200;

	@Parameter(names={"--score"}, description="What type of score system to use for peaks. Can be used for ranking peaks. " +
									 		  "Options: <max || ave || med || fc || zscore || all>  Default = max")
	public String scoreSys = "max";

	@Parameter(names={"--bgscore"}, arity=1, description="Whether to add the HMMR score to each state annotation in bedgraph. Note: this adds considerable time. Default = False")
	public boolean bgScore = false;

	@Parameter(names={"-k", "--kmeans"}, description="Number of States in the model. Default = 3. If not k=3, recommend NOT calling peaks, use bedgraph")
	public int k = 3;

	@Parameter(names={"--trim"}, description="How many signals from the end to trim off (i.e. starting with tri signal then di etc). " +
											 "This may be useful if your data doesn't contain many large fragments. Default = 0")
	public int trim = 0;

	@Parameter(names={"-t", "--training"}, converter=FileConverter.class, description="BED file of training regions to use for training model, instead of foldchange settings")
	public File trainingRegions = null;

	@Parameter(names={"--window"}, description="Size of the bins to split the genome into for Viterbi decoding." +
											   "To save memory, the genome is split into <int> long bins and viterbi decoding occurs across each bin. " +
											   "Default = 25000000. Note: For machines with limited memory, it is recommended to reduce the size of the bins.")
	public int vitWindow = 25000000;

	@Parameter(names={"--model"}, converter=FileConverter.class, description="Binary model file (generated from previous HMMR run) to use instead of creating new one")
	public File modelFile = null;

	@Parameter(names={"--modelonly"}, arity=1, description="Whether or not to stop the program after generating model. Default = false")
	public boolean stopAfterModel = false;

	@Parameter(names={"--printTracks"}, arity=1, description="Whether or not to print the decomposed HMMRATAC signal tracks. " +
															 "Tracks will be labeled as Output_NFR.bedgraph, Output_Mono.bedgraph etc. Default = false")
	public boolean printHMMRTracks = false;

	@Parameter(names={"--maxTrain"}, description="Maximum number of training regions to use. Default = 1000")
	public int maxTrain = 1000;

	@Parameter(names={"--removeDuplicates"}, arity=1, description="Whether or not to remove duplicate reads from analysis. Default = true")
	public boolean rmDup = true;

	@Parameter(names={"--printExclude"}, arity=1, description="Whether to output excluded regions into Output_exclude.bed. Default = false")
	public boolean printExclude = false;

	@Parameter(names={"--printTrain"}, arity=1, description="Whether to output training regions into Output_training.bed. Default = true")
	public boolean printTrain = true;

	@Parameter(names={"--randomSeed"}, description="Seed to set for random sampling of training regions. Default is 10151")
	public long randomTrainSeed = 10151;

	@Parameter(names={"--threshold"}, description="threshold for reporting peaks. Only peaks who's score is >= this value will be reported.")
	public double threshold = 30;

	@Parameter(names={"-h", "--help"}, help=true)
	private boolean help;
}
