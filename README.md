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
The source files are uploaded for users to see how the program works.  They do not contain the required manifest file needed to actually
run.  For details on HOW to run HMMRATAC, see HMMRATAC_Guide.txt, which contains a thorough runthrough of all parameters, output files and input
requirements and troubleshooting.
