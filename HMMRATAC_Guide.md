# HMMRATAC Guide

Running the program:

```java -jar HMMRATAC_V1.2_exe.jar -b <SortedBAM> -i <BAMIndex> -g <GenomeStatsFile> <options>```

Version number may change.

Please download the executable file from
https://github.com/LiuLabUB/HMMRATAC/releases to run the
program. Source files do not contain manifest file needed to run.

### Outline:
1. [Creating required files](#creating-required-files)
2. [Commandline options](#commandline-options)
3. [Output files](#output-files)
4. [Troubleshooting](#troubleshooting)

## Creating required files

1. Sorted BAM file containing the alignments of ATAC-seq reads

   The first file that is necessary is a sorted BAM file containing
   the pair-end ATAC-seq reads.  Typically, an alignment algorithm
   (such as bwa/bowtie2/minimap2) will output a SAM file as default.
   The first step is to convert this SAM file into a BAM file.  This
   is best accomplished using ```samtools```.  Samtools requires a SAM
   file to convert into a BAM file. Use the ```samtools view```
   command as follows, assuming the SAM file is named as
   ```atac.sam```:

   ```samtools view –bS -o atac.bam atac.sam```

   **Note** the SAM file should have header lines. If not, check the
     step3 below to get the genome information file ready
     ```genome.info```(chromosome lengths) first, and run
     ```samtools``` as:

   ```samtools view -bS -t genome.info -o atac.bam atac.sam```
   
   The next step is to sort the BAM file. This can also be
   accomplished using samtools.  Use the ```samtools sort``` command
   as follows (example using 4 threads to sort, tweak it as
   neccessary):

   ```samtools sort -@ 4 -o atac.sorted.bam atac.bam```

   It is possible to remove duplicates or low Mapping quality reads
   before inputting into HMMRATAC, although HMMRATAC will perform
   these functions by default.  By default HMMRATAC will remove
   duplicate reads (and this function is currently hard-coded and
   cannot be turned off).  HMMRATAC also removes those reads whose
   MapQ (mapping quality scores) are below 30.  This value is
   changeable through command line option (using the ```–q``` or
   ```--minmapq``` option)

   If you have multiple replicates and you choose to merge replicates
   prior to running HMMRATAC, this can also be accomplished using
   samtools. Use the ```samtools merge``` command as follows:

   ```samtools merge –n atac_merged.bam atac_rep1.sorted.bam atac_rep2.sorted.bam ... atac_repN.sorted.bam```

   If you merge BAM files sorted by the same ```samtools``` tool, the
   merged file ```atac_merged.bam``` should *be sorted already*. If
   uncertain, you can always re-sort the resulting merged file before
   proceeding.

   ```samtools sort -@ 4 -o atac_merged.sorted.bam atac_merged.bam```

   Now, either rename ```atac.sorted.bam``` if you have only one
   replicate, or ```atac_merged.bam``` if you have multiple
   replicates, or ```atac_merged.sorted.bam``` if an extra sorting is
   performed on merged file, to ```atac.forHMMRATAC.bam```.

2. The index file for the BAM.

   The next required file that HMMRATAC needs is an index file for the
   sorted BAM file.  This file can also be easily created with
   samtools.  Use the ```samtools index``` command as follows:

   ``` samtools index atac.forHMMRATAC.bam```

   You will have a file named ```atac.forHMMRATAC.bam.bai``` in the
   same working directory.

3. Genome information file for chromosome sizes.

   The third required file is a genome information file containing
   chromosome sizes.  This file can be downloaded from numerous online
   sources, including UCSC genome browser’s website.  This file is a
   two column, tab-delimited file containing chromosome name and size
   in the following format: <chromName><TAB><chromSIZE>

** Please Note: The genome file must NOT contain a header line

   It is possible to retrieve this file using UCSC Genome Browser’s
   MySQL database.  To retrieve the file, for H. sapiens, use the
   following command:

   ```mysql –-user=genome --host=genome-mysql.cse.ucsc.edu –A –e \ “select chrom, size from hg19.chromInfo” > genome.info ```

   Alternatively, if your SAM file has header lines and you have
   successfully generated a BAM file following step 1 and 2 without
   input of genome information. You can extract genome information
   from BAM as:

   ```samtools view -H atac.forHMMRATAC.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info```

## Commandline Options

### Required Parameters: 

* ```-b , --bam <BAM>```

   This is the sorted BAM file created as described in
   section 1.1.

* ```-i , --index <BAI>```

   This is the index file for the sorted BAM file
   created as described in [section 1](#creating-required-files).

* ```-g , --genome <GenomeFile>```

   This is the genome file created or
   downloaded as described in [section 1](#creating-required-files).

### Optional Parameters:

* ```-m , --means <Double>```

   Comma separated list of initial mean values for the fragment
   distribution, used to create the signal tracks.  The default values
   are 50,200,400,600.  It is recommended to change these values if
   you are working with non-human species, as other species have
   different nucleosome spacing.  These values will be updated with EM
   training, unless the ```-f``` option is set to false.  Also note that the
   first value, corresponding to the short read distribution, is NOT
   updated with EM training and will remain as it is set with this
   option.  It is recommended to set this first value to the read
   length, if the read length is under 100bp.

* ```-s , --stddev <Double>```

   Comma separated list of initial standard deviations values for the
   fragment distributions., used to create the signal tracks.  The
   default values are 20,20,20,20.  These values may also need to be
   updated for non-human species, although it may not be necessary, as
   the EM is robust in handling variances.  These values are also
   updated with EM training.  Note: the first value is not updated, as
   the short distribution uses an Exponential distribution that only
   uses the mean value, so the first value is meaningless.

* ```-f , --fragem <True || False>```

   Boolean to determine whether fragment EM training is to occur or
   not.  The default is true.  If this option is set to false, the
   initial values described with ```-m``` and ```-s``` are used as the
   parameters of the mixture model.  This is generally not
   recommended.  Setting this value to false can decrease the total
   runtime, but may result in a less accurate model.  If the data has
   been run already, and it is desired to re-run it using different
   reporting thresholds, you could set this option to false, provided
   that you reset the initial parameters to the updated ones created
   by the previous run.  These values are recorded in the .log file
   described later.

* ```-q , --minmapq <int>```

   This is the minimum mapping quality score for reads to be used in
   creating the signal tracks.  The default is 30.  Although this can
   be set to higher or lower values, it is generally not recommended.
   If few meet this threshold, that would indicate that the assay
   itself was compromised, and setting a lower value would likely
   still result in errors.

* ```-u , --upper <int>```

   Upper limit on fold change range for choosing training regions.
   HMMRATAC chooses training regions by finding genomic loci whose fold
   change above genomic background is within a certain range.  This
   option sets the upper limit of this range.  Default is 20.
   Generally speaking, the higher this range is, the more stringent
   the resulting model, and the lower the range, the less stringent
   the model.  If recall (sensitivity) is important for you, it is
   recommended to set lower ranges, while if precision or accuracy is
   more important, it is recommended to set higher ranges.  Once the
   regions within a certain range are found, HMMRATAC extends the regions
   by +/- 5KB and uses these extended regions to train the model.  The
   search ends once a maximum of 1000 regions are identified.

* ```-l , --lower <int>```

   Lower limit on fold change range for choosing
   training regions.  HMMRATAC chooses training regions by finding genomic
   loci whose fold change above genomic background is within a certain
   range.  This option sets the lower limit of this range.  Default is
   10.  Generally speaking, the higher this range is, the more
   stringent the resulting model, and the lower the range, the less
   stringent the model.  If recall (sensitivity) is important for you,
   it is recommended to set lower ranges, while if precision or
   accuracy is more important, it is recommended to set higher ranges.
   Once the regions within a certain range are found, HMMRATAC extends
   regions by +- 5KB and uses these extended regions to train the
   model.  The search ends once a maximum of 1000 regions are
   identified.

* ```-z , --zscore <int>```

   Zscored read depth to mask during Viterbi decoding.  Default is
   100.  In our experience, Viterbi has trouble decoding regions that
   Have very high read coverage.  If encountered, Viterbi will either
   call a Empty block, where no state is called, or will call one
   large open region.  To avoid this problem, HMMRATAC will skip over any
   regions whose centered Read coverage is equal to or greater than
   this value.  If the report peaks Option is set (-p True – described
   below), these regions are added back to The peak file and are
   labeled as “High_Coverage_Peak_#”.  These regions Can therefore be
   examined further, if desired.

* ```-o , --output <Name>```

   The prefix for all output files.  Default is “NA”.

* ```-e , --blacklist <BED>```

   BED file of regions to exclude from model creation and genome
   decoding.  These could be previously annotated blacklist regions.
   It is always recommended to include a blacklisted region list.

* ```-p , --peaks <True || False>```

   Boolean to determine whether or not to report a peak file in BED
   format.  Default is True.  The resulting file will be in gappedPeak
   format and will be described in detail below.

* ```-k , --kmeans <int>```

   Number of states in the model.  Default is 3.  It is generally not
   recommended to change this setting.  All of our tests and
   comparisons were done using the default setting.  It may be
   reasonable to set to a lower number of states when using 500 cell
   per-replicate data, but this is untested.  If this setting is
   changed, it is generally recommended to not report peaks, but only
   to report the genome-wide state annotation file.  This can then be
   interpreted manually.

* ```-t , --training <BED>```

   BED file of training regions to use instead of using fold-change
   ranges.  If this option is used, it is recommended to only use 1000
   regions and to extend them by +/- 5KB, as HMMRATAC would do for
   fold-change identified training regions.

* ```--bedgraph <True || False>```

   Boolean to determine whether a genome-wide state annotation
   bedgraph file should be reported.  Default is False.

* ```--minlen <int>```

   Minimum length of an open region state to call a peak.  Note: that
   the ```-p``` option must be set to true. Default is 200.

* ```--score <max || ave || med || fc || zscore || all>```

   What type of score system to use for scoring peaks.  “Max” refers
   to the maximum number of reads mapping to the open region.  “Med”
   refers to the median number of reads mapping to the open region.
   “Ave” refers to the average number of reads mapping to the open
   region.  “FC” refers to fold-change and is the average number of
   reads mapping to the open region divided by the genome average.
   “Zscore” refers to the average number of reads mapping to the open
   region minus the genome average divided by the genomic standard
   deviation.  “All” reports all these scores separated by a “_”
   (underscore). The order for "all" is MAX_MEAN_MEDIAN_ZSCORE_FC
   Default is “max”.

* ```--bgscore <True || False>```

   Boolean to determine whether to add a score to every state
   annotation in a bedgraph file.  Note that ```--bedgraph``` has to be set
   true.  This adds considerable time to the program and is generally
   not recommended. Default is false.

* ```--trim <int>```

   How many signals, or distributions, to trim from the signal tracks.
   It trims from the end of the matrix (IE 1 means trim the tri signal
   track, 2 means trim the tri  and di signal tracks…etc).  This could
   be useful if your data  doesn’t contain many longer fragments, such
   as  size selected  data.  Recommendations:  fragments <=  500bp set
   ```--trim 1```; fragments <= 250bp set ```--trim 2```; fragments <=
   150 set ```--trim 3```.  Default is 0.

* ```--window <int>```

   Size of the bins to split the genome into for Viterbi decoding.  To
   save memory, HMMRATAC splits the genome into these sized bins and
   Performs Viterbi decoding on these bins separately.  Default is
   25000000.  It may be necessary to reduce the size of these bins, if
   running HMMRATAC on A desktop or a machine with limited memory.  Most
   desktops should handle A bin size of 1/10 the default.  A java heap
   space runtime error, is common When the bin size is too large for
   the machine (see [section 4](#troubleshooting)).

* ```--model <File>```

   This model (binary model generated by previous HMMRATAC run, suffixed
   with .model) will be used to decode the genome rather than building
   a new model using training regions (either provided by user with
   the ```-t``` option or created by HMMRATAC using the ```-u``` and ```-l```
   options).  This can be used in conjunction with the ```--modelonly```
   option (create the model, inspect it and run HMMRATAC with that model).

* ```--modelonly <True || False>```

   Boolean to determine if HMMRATAC should quit after generating the
   model.  This is a helpful option when trying to determine the best
   parameters for creating the model.  Default = false.

* ```--maxTrain <int>```

   Maximum number of training regions to use during model building.
   It is possible that the default number of training regions (1000) is 
   too many and will cause memory issues for smaller machines. Setting 
   this parameter to 1/2 of the default (ie. 500) could help with this 
   issue.

* ```-h , --help``` This flag prints a help message and exits the program.

## Output Files
	
1. Name.log

  This is a log file generated as HMMRATAC runs.  It contains all inputted
  arguments, the results of the fragment distribution EM algorithm,
  and various status updates as HMMRATAC progresses, as well as a textual
  description of the generated model.  When finished, it also reports
  the total amount of actual (real-world) time it took for HMMRATAC to
  run.

2. Name.bedgraph

  This is the genome-wide state annotation file created if the
  ```--bedgraph``` option is set to true.  It is a four column,
  tab-delimited file, with the following fields:

  1. Chromosome Name

  2. Region Start (zero-based)

  3. Region Stop (zero-based)

  4. State annotation. This represents what state the particular
     region is assigned to, corresponding to the created model.  The
     state name is prefixed with the letter “E”, to make extractions
     easier.  For instance, to retrieve all regions that were
     classified as state 1, type:

     ```grep E1 Name.bedgraph > State1_regions.bed```

3. Name.model

  This file contains the model that was generated by HMMRATAC and was
  used to decode the genome.  It is always recommended to check this
  model after running HMMRATAC.  Please note that this file is a
  binary file that can only be read by HMMRATAC.  To see the actual
  textual description of the model, see the log file.  If HMMRATAC
  runs to completion but produces poor results, it is usually the
  result of a faulty model.  For instance, if the model shows “NaN”
  for all of the parameters, such as transition or emission
  probabilities, this usually means that something went wrong with
  either the choice of training regions or generating the signal
  tracks.  This can occur with too low or too high of a fold-change
  range for choosing training regions, if the BED file of training
  regions was faulty or if certain signals need to be trimmed off (due
  to small numbers of larger fragments).  Generally speaking, this
  poor model occurs when there is not enough separation in the signal
  tracks between the states.  Common fixes include, trimming the
  signal tracks (option ```--trim```), changing the number of states
  (option ```-k```), choosing a different fold change range for
  training site identification (option ```-u``` and ```-l```) or
  changing the list of previously annotated training regions (option
  ```-t```).  Additionally, checking the model ensures that the
  predicted model is created.  Although rare, it is possible for a
  different state to be better indicative of the open state, than is
  normally the case (HMMRATAC sorts the model so the last or highest
  numbered state should be the open state).  If this happens, then the
  resulting peak file may be incorrect.  In such cases, it is
  recommended to report a bedgraph file and extract the proper state
  (of certain lengths and with the accompanying flanking regions) and
  manually create a peak BED file.  We have tested numerous ATAC-seq
  datasets, using our default parameters and we generally don’t
  encounter problems, unless we radically change those settings.

4. Name_peaks.gappedPeak

  This is the standard peak file reported by HMMRATAC by default (when
  the ```-p , --peaks``` option is set to true).  The first line in
  the file is a track line for genome browsers.  After that, it is a
  15 column, tab-delimited file containg the following fields:

  1. Chromosome Name

  2. Peak Start (zero-based). This is the beginning of the regulatory
region that HMMRATAC identifies as a peak.  This includes the flanking
nucleosome regions.  Therefore, this position represents the start of
the upstream nucleosome.  If there is no upstream nucleosome, this
position represents the start of the open state.

  3. Peak Stop (zero-based). This is the end of the regulatory region
that HMMRATAC identifies as a peak.  This includes the flanking
nucleosome regions.  Therefore, this position represents the end of
the downstream nucleosome.  If there is no downstream nucleosome, this
position represents the end of the open state.

  4. Peak Name. Unique name for the peak, in the format “Peak_#”.
For the high coverage regions that are excluded from Viterbi decoding
(described in [Section 2](#commandline-options), under the ```-z ,
--zscore``` option), the name is in the format “HighCoveragePeak_#”.
This allows the high coverage peaks to be easily identified and
extracted.

  5. Unused, denoted with “.”

  6. Unused, denoted  with “.”

  7. Open state start (zero-based).  This is the genomic position
where the open state begins.  It is therefore possible to use only the
open regions, rather than the regulatory regions that HMMRATAC reports
by default.

  8. Open state stop (zero-based).  This is the genomic position where
the open state ends.  It is therefore possible to use only the open
regions, rather than the regulatory regions that HMMRATAC reports by
default.

  9. Color code for display. Set to ```255,0,0```.  This would create a
dark-blue color on a genome browser.

  10. The number of sub-regions in the peak region.  For standard
peaks this is set to 3, indicating the open state region and the two
flanking nucleosome state regions.  For high coverage peaks, this is
set to 1.  Peaks that don’t have upstream and/or downstream
nucleosomes will have different values as well.

  11. Comma separated list of sub-region lengths.  The first and last
values are always set to one, making visualization easier.  The middle
value is the length of the open state region.  For high coverage
peaks, this is only a single value, denoting the length of the total
region.  If the peak lacks upstream and/or downstream nucleosomes,
these values will reflect that.

  12. Comma separated list of the (zero-based) starts of each
sub-region, relative to the total region’s start.  For instance, the
first value is always 0, meaning 0 bp from the region start (column
2).  The second value is the distance from the region start (column 2)
to the open region start (column 7), etc.

  13. Peak Score. This is the score for the peak, as determined by the
  scoring system option described in [section
  2](#commandline-options) (option ```--score```).

  14. Unused, denoted as -1.

  15. Unused, denomted as -1.

5. Name_summits.bed

  This is a BED file containing the summits of each HMMRATAC peak.
  Summits are determined by finding the position within the open state
  region, whose Gaussian-smoothed read coverage is the maximum over
  the entire region.  This is only outputted if a peak file is also
  outputted.  If motif analysis is being conducted, it is recommended
  to use these summits as the points of interest, as the summits tend
  to be better indicative of TF binding as opposed to peak centers.
  This is a four column, tab-delimited BED file containing the
  following fields:

  1. Chromosome Name

  2. Summit Start (zero-based)

  3. Summit Stop (zero-based)

  4. Summit Name.  This is the same name as the corresponding Peak
  name used in the peak file (column 4 of the Name_peaks.gappedPeak
  file).

  5. Peak Score. Same score as in the .gappedPeak file.

## Troubleshooting

1. Java heap space error:

   The most likely cause of this error is that the genome "blocks" are
   too large. In order to save memory, HMMRATAC splits the genome into
   "blocks" and then annotates each position within the block as one
   of the states of the model. The default size of these blocks is
   25MB. We have found, that on personal computers, this may be too
   large to handle, resulting in the heap space error. The best fix to
   this issue is to reduce the size of the block, using the
   ```--window``` option (See [section 2](#commandline-options)).
   It is recommended to reduce this value to 1/10 of the default, ie
   2.5MB, for personal computers.  If the error continues, you can
   lower it further.

2. HMMRATAC produces very large peaks (> 25KB) or results in large
gaps with no peaks:

   This isssue is related to the Viterbi algorithm.  In our
   experience, Viterbi has difficulty when encountering very high
   coverage regions. When it does, one of two things can
   happen. Either viterbi will call every position, from the high
   coverage region to the end of the block (block is explained in
   "Java heap space error" section) as a peak (thereby resulting in
   extremely large peaks) or as some other state (resulting in large
   gaps with no called peaks).  HMMRATAC corrects for this behavior by
   skipping over high coverage regions, allowing viterbi to continue
   correctly.  HMMRATAC determines what is or isn't a high coverage
   region based on the ```-z``` option (see [section
   2](#commandline-options)).  Lowering the ```-z``` option, will
   make HMMRATAC skip more high coverage regions (by lowering the
   threshold to identify such regions). It should be noted that all
   skipped high coverage regions are added back to the gappedPeak
   output file and are annotated as "High Coverge Peaks".

3. HMMRATAC produces very few peaks (< 1000) or doesn't identify any
peaks or the model shows NaN for all values:

   This is almost always due to a faulty model.  Check the log file
   (see [section 3](#output-files)) and see if the model is valid. A non-valid
   model will have NaN as the values of the initial, transition and
   emission probabilities.  This typically occurs when there is not
   enough separation in the values between states. As an example
   (using a simplified model), say you had a model with 3 states and a
   single observation value.  The average value for state 1 is 1, the
   average for state 2 is 2 and the average for state 3 is 3.  If a
   value of 2.5 occurred, then there is an equal chance that the value
   is assigned to state 2 or state 3.  When this happens, HMMRATAC (or
   really any model based algorithm) won’t be able to decide which
   state to assign it to.  Now we remake the model with a more
   conservative training set and the values become 1, 5 and 10.  Now
   2.5 clearly can be assigned to state 1. Therefore, an easy way to
   correct this is to use more stringent settings in creating the
   model. This is accomplished with the ```-l``` and ```-u``` options
   (see [section 2](#commandline-options)). The default setting for
   ```-u``` is 20 and for ```-l``` is 10. A common fix is to increase
   the ```-u``` to 30 or 40. This will allow stronger peaks into the
   training set, thereby changing the emission values of each state
   and correcting the problem (not enough separation between states).
   This is the most common error with HMMRATAC.

4. The dataset should not be empty ERROR. This error has been
encountered on occasion, and is almost always the result of user
error.

   This is saying that the data matrix created by HMMRATAC, containing
   the values for nucleosome free and nucleosome enriched fragments,
   was empty. If you set the upper and lower fold change bounds for
   choosing training regions (```-u``` and ```-l``` option) to values
   that prevent ANY training sites from being identified (i.e. ```-u
   10 -l 20```), or if you provide an empty BED file of training
   regions, then the data matrix will be empty and this error will
   occur. Alternatively, if the BAM files use one type of chromsome
   annotation (ie. chr1) while the genome file uses a different
   chromosome annotation (i.e. 1) or vice versa, then this discrepancy
   will also result in an empty data matrix. Finally, it may be that
   there were no identified training regions even with reasonable
   settings, than it may be neccessary to choose a different method to
   choose training sites, i.e. to set different ```-u``` and ```-l```
   settings.

5. Tips for checking model.

   The creation of a proper model is critial for HMMRATAC.  It is
   always recommended to check the log file and visually examine the
   model to confirm that a proper one has been created.  The open
   state should have the highest emission parameters for all 4 signal
   components, the nucleosome state should have the second highest
   emission parameters for all 4 signal components and the background
   state should have the lowest emission parameters for all 4 signal
   components.  If this is not the case, the model is likely bad.  An
   example is a case where the open state has the highest NFR and
   mononucleosome emission values, but the second highest dinucleosome
   and trinucleosome emission values.  This often happens when there
   isn't enough separation in the values between states (similar to
   item 3 above).  A typcial fix would be to re-run the software with
   a more stringent set of parameters (```-u``` and ```-l```, similar
   to item 3 above).  Another possible issue is that the log file
   shows that the kmeans model is identical to the final model.  In
   the event that Baum Welch training fails, HMMRATAC will use the
   original Kmeans model in order to produce some result.  Because the
   kmeans model doesn't include transition parameters, this is also
   not ideal. The way to solve this issue, is the same as above,
   namely setting mote stringent parameters to create more separation
   in the values between states.
   
6. HMMRATAC produces too many peaks. What to keep?

   By default, HMMRATAC will report all peak regions that match the learned 
   model. This includes weaker peaks that have the characteristic patterns of
   nucleosome-free and nucleosome-enriched fragments but whose read depth is low.
   It is often necessary to filter the resulting peaks by their score, in order 
   to eliminate these weak peaks, if that is desired. As an example, say you called 
   peaks and wanted to keep only those whose score (defined by the --score option - 
   default is maximum reads in the peak) is greater than or equal to 10. You 
   could use the following command to filter the peak file:
   
    ``` awk -v OFS="\t" '$13>=10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak```
   
   To filter the summit file by the same threshold:
   
   ```awk -v OFS="\t" '$5>=10 {print}' NAME_summits.bed > NAME.filteredSummits.bed```
   
   Generally speaking, the higher the threshold that is used for filtering, the higher
   the precision will be and the lower the sensitivity will be.  A lower score will
   result in higher sensitivity and lower precision.
   
