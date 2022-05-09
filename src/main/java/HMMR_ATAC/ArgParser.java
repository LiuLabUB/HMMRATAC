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

import com.beust.jcommander.*;

import java.io.File;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;


public class ArgParser {
	
	private static class FileConverter implements IStringConverter<File> {
		@Override
		public File convert(String value) {
			return new File(value);
		}
	}
	
	public static class MyUsageFormatter extends DefaultUsageFormatter {
		public MyUsageFormatter(JCommander commander) {
			super(commander);
		}
		
		private static String newLineAndIndent(int indent) {
			return "\n" + s(indent);
		}
		
		public void appendAllParametersDetails(StringBuilder out, int indentCount, String indent,
											   List<ParameterDescription> sortedParameters) {
			if (sortedParameters.size() > 0) {
				out.append(indent).append("  Options:\n");
			}
			
			for (ParameterDescription pd : sortedParameters) {
				WrappedParameter parameter = pd.getParameter();
				String description = pd.getDescription();
				boolean hasDescription = !description.isEmpty();
				
				// First line, command name
				out.append(indent)
						.append("  ")
						.append(parameter.required() ? "* " : "  ")
						.append(pd.getNames())
						.append("\n");
				
				if (hasDescription) {
					wrapDescription(out, indentCount, s(indentCount) + description);
				}
				Object def = pd.getDefault();
				
				if (pd.isDynamicParameter()) {
					String syntax = "Syntax: " + parameter.names()[0] + "key" + parameter.getAssignment() + "value";
					
					if (hasDescription) {
						out.append(newLineAndIndent(indentCount));
					} else {
						out.append(s(indentCount));
					}
					out.append(syntax);
				}
				
				if (def != null && !pd.isHelp()) {
					String displayedDef;
					if (def instanceof double[]) {
						displayedDef = Arrays.toString((double[]) def);
					} else if (Strings.isStringEmpty(def.toString())) {
						displayedDef = "<empty string>";
					} else {
						displayedDef = def.toString();
					}
					String defaultText = "Default: " + (parameter.password() ? "********" : displayedDef);
					
					if (hasDescription) {
						out.append(newLineAndIndent(indentCount));
					} else {
						out.append(s(indentCount));
					}
					out.append(defaultText);
				}
				Class<?> type = pd.getParameterized().getType();
				
				if (type.isEnum()) {
					String valueList = EnumSet.allOf((Class<? extends Enum>) type).toString();
					String possibleValues = "Possible Values: " + valueList;
					
					// Prevent duplicate values list, since it is set as 'Options: [values]' if the description
					// of an enum field is empty in ParameterDescription#init(..)
					if (!description.contains("Options: " + valueList)) {
						if (hasDescription) {
							out.append(newLineAndIndent(indentCount));
						} else {
							out.append(s(indentCount));
						}
						out.append(possibleValues);
					}
				}
				out.append("\n");
			}
		}
	}
	
	//Required inputs
	
	@Parameter(names = {"-b", "--bam"}, converter = FileConverter.class, required = true, order = 0,
			description = "Sorted BAM file containing the ATAC-seq reads")
	public File bam;
	
	@Parameter(names = {"-i", "--index"}, converter = FileConverter.class, required = true, order = 0,
			description = "Index file for the sorted BAM File")
	public File index;
	
	@Parameter(names = {"-g", "--genome"}, converter = FileConverter.class, required = true, order = 0,
			description = "Two column, tab delimited file containing genome size stats")
	public File genomeFile;
	
	//Optional Inputs
	
	@Parameter(names = {"-m", "--means"}, arity = 4, order = 1,
			description = "Comma separated list of initial mean values for the fragment distribution.")
	public double[] means = {50.0, 200.0, 400.0, 600.0};
	
	@Parameter(names = {"-s", "--stddev"}, arity = 4, order = 1,
			description = "Comma separated list of initial standard deviation values for fragment distribution.")
	public double[] stddevs = {20.0, 20.0, 20.0, 20.0};
	
	@Parameter(names = {"-f", "--fragem"}, arity = 1, order = 1,
			description = "Whether to perform EM training on the fragment distribution.")
	public boolean fragEM = true;
	
	@Parameter(names = {"-q", "--minmapq"}, order = 1,
			description = "Minimum mapping quality of reads to keep.")
	public int minMapQ = 30;
	
	@Parameter(names = {"-u", "--upper"}, order = 1,
			description = "Upper limit on fold change range for choosing training sites.")
	public int upper = 20;
	
	@Parameter(names = {"-l", "--lower"}, order = 1,
			description = "Lower limit on fold change range for choosing training sites.")
	public int lower = 10;
	
	@Parameter(names = {"-z", "--zscore"}, order = 1,
			description = "Zscored read depth to mask during Viterbi decoding.")
	public int zscore = 100;
	
	@Parameter(names = {"-o", "--output"}, order = 1,
			description = "Name for output files.")
	public String output = "NA";
	
	@Parameter(names = {"-p", "--peaks"}, arity = 1, order = 1,
			description = "Whether to report peaks in bed format.")
	public boolean peaks = true;
	
	@Parameter(names = {"--bedgraph"}, arity = 1, order = 2,
			description = "Whether to report whole genome bedgraph of all state annotations.")
	public boolean bg = false;
	
	@Parameter(names = {"-e", "--blacklist"}, converter = FileConverter.class, order = 1,
			description = "bed file of blacklisted regions to exclude")
	public File blacklist = null;
	
	@Parameter(names = {"--minlen"}, order = 2,
			description = "Minimum length of open region to call peak. Note: -p , --peaks must be set.")
	public int minLength = 200;
	
	@Parameter(names = {"--score"}, order = 2,
			description = "What type of score system to use for peaks. Can be used for ranking peaks. " +
					"Options: <max || ave || med || fc || zscore || all> ")
	public String scoreSys = "max";
	
	@Parameter(names = {"--bgscore"}, arity = 1, order = 2,
			description = "Whether to add the HMMR score to each state annotation in bedgraph. Note: this adds considerable time.")
	public boolean bgScore = false;
	
	@Parameter(names = {"-k", "--kmeans"}, order = 1,
			description = "Number of States in the model. If not k=3, recommend NOT calling peaks, use bedgraph")
	public int k = 3;
	
	@Parameter(names = {"--trim"}, order = 2,
			description = "How many signals from the end to trim off (i.e. starting with tri signal then di etc). " +
					"This may be useful if your data doesn't contain many large fragments.")
	public int trim = 0;
	
	@Parameter(names = {"-t", "--training"}, converter = FileConverter.class, order = 1,
			description = "BED file of training regions to use for training model, instead of foldchange settings")
	public File trainingRegions = null;
	
	@Parameter(names = {"--window"}, order = 2,
			description = "Size of the bins to split the genome into for Viterbi decoding." +
					"To save memory, the genome is split into <int> long bins and viterbi decoding occurs across each bin. " +
					"Note: For machines with limited memory, it is recommended to reduce the size of the bins.")
	public int vitWindow = 25000000;
	
	@Parameter(names = {"--model"}, converter = FileConverter.class, order = 2,
			description = "Binary model file (generated from previous HMMR run) to use instead of creating new one")
	public File modelFile = null;
	
	@Parameter(names = {"--modelonly"}, arity = 1, order = 2,
			description = "Whether or not to stop the program after generating model.")
	public boolean stopAfterModel = false;
	
	@Parameter(names = {"--printTracks"}, arity = 1, order = 2,
			description = "Whether or not to print the decomposed HMMRATAC signal tracks. " +
					"Tracks will be labeled as Output_NFR.bedgraph, Output_Mono.bedgraph etc.")
	public boolean printHMMRTracks = false;
	
	@Parameter(names = {"--maxTrain"}, order = 2,
			description = "Maximum number of training regions to use.")
	public int maxTrain = 1000;
	
	@Parameter(names = {"--removeDuplicates"}, arity = 1, order = 2,
			description = "Whether or not to remove duplicate reads from analysis.")
	public boolean rmDup = true;
	
	@Parameter(names = {"--printExclude"}, arity = 1, order = 2,
			description = "Whether to output excluded regions into Output_exclude.bed.")
	public boolean printExclude = false;
	
	@Parameter(names = {"--printTrain"}, arity = 1, order = 2,
			description = "Whether to output training regions into Output_training.bed.")
	public boolean printTrain = true;
	
	@Parameter(names = {"--randomSeed"}, order = 2,
			description = "Seed to set for random sampling of training regions.")
	public long randomTrainSeed = 10151;
	
	@Parameter(names = {"--threshold"}, order = 2,
			description = "threshold for reporting peaks. " +
					"Only peaks whose score is greater than or equal to the threshold will be reported.")
	public double threshold = 30;
	
	@Parameter(names = {"-h", "--help"}, help = true,
			description = "Print this help message and exit.")
	public boolean help = false;
}
