# HMMRATAC

**Quick Start**

Assume that you have a BAM file from aligner such as ```bwa mem``` named ```ATACseq.bam```.

1. Sort the BAM file to get a ```ATACseq.sorted.bam``` file:

   ```samtools sort ATACseq.bam  -o ATACseq.sorted.bam```

2. Make index from the BAM file to get a ```ATACseq.sorted.bam.bai``` file:

   ```samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai```

3. Make genome information (chromosome sizes) from the BAM file to get a ```genome.info``` file:

   ```samtools view -H ATACseq.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info ```

4. Run HMMRATAC on the sorted BAM ```ATACseq.sorted.bam```, the BAM index file ```ATACseq.sorted.bam.bai```, and the genome information file ```genome.info```:

   ```java -jar HMMRATAC_V1.2.4_exe.jar -b ATACseq.sorted.bam -i ATACseq.sorted.bam.bai -g genome.info```

Samtools can be downloaded here: http://www.htslib.org/download/

Be sure to run HMMRATAC using the executable file, found here: 
https://github.com/LiuLabUB/HMMRATAC/releases
For details on HOW to run HMMRATAC, see HMMRATAC_Guide.md, which contains a thorough runthrough of all parameters, output files and input
requirements and troubleshooting.

HMMRATAC requires paired-end data. Single-end data will not work. HMMRATAC is designed to process ATAC-seq data that hasn't undergone
any size selection, either physical or in silico. This should be standard practice for any ATAC-seq analysis. Size selected data can be
processed by HMMRATAC (see [HMMRATAC_Guide.md](./HMMRATAC_Guide.md#commandline-options) on ```--trim``` option). 

If you use HMMRATAC in your research, please cite the following paper:

Evan D. Tarbell and Tao Liu, [HMMRATAC: a Hidden Markov ModeleR for ATAC-seq](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz533/5519166), Nucleic Acids Res. 2019 Jun 14. doi: 10.1093/nar/gkz533
