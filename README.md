# HMMRATAC

Quick Start:

$ samtools sort ExampleFile.bam  -o ExampleFile.sorted.bam

$ samtools index ExampleFile.sorted.bam ExampleFile.sorted.bam.bai

$ samtools view -H ExampleFile.sorted.bam | grep SQ | cut -f 2-3 | cut -d ':' -f 2,3 | cut -d 'L' -f 1 > tmp

$ samtools view -H ExampleFile.sorted.bam | grep SQ | cut -f 2-3 | cut -d ':' -f 2,3 | cut -d ':' -f 2 > tmp2

$ paste tmp tmp2 > hg.genome; rm tmp tmp2

$ java -jar HMMRATAC_V1.2_exe.jar -b ExampleFile.sorted.bam -i ExampleFile.sorted.bam.bai -g hg.genome

**NOTE: Earlier versions of HMMRATAC require a bigwig file. See HMMRATAC_Guide.txt for more detail**

Samtools can be downloaded here: http://www.htslib.org/download/

Be sure to run HMMRATAC using the executable file, found here: 
https://github.com/LiuLabUB/HMMRATAC/releases
For details on HOW to run HMMRATAC, see HMMRATAC_Guide.txt, which contains a thorough runthrough of all parameters, output files and input
requirements and troubleshooting.

HMMRATAC requires paired-end data. Single-end data will not work. HMMRATAC is designed to process ATAC-seq data that hasn't undergone
any size selection.

If you use HMMRATAC in your research, please cite the following paper:

Evan D. Tarbell and Tao Liu, "HMMRATAC, a Hidden Markov ModeleR for ATAC-seq", bioRxiv 306621, Pages 1-24, 2018. https://doi.org/10.1101/306621 
