import sys
import commands as cm

if len(sys.argv) < 4:
	print "Usage: python Make_HMMR_Files.py <UnsortedBAMFile> <GenomeStatsFile> <Prefix>"
	sys.exit()
#Absolute path to required tools
#Please provide the correct path IE:
#	samtools = "/path/to/samtools"
#	Note: quotation marks are required

samtools = "samtools"
bedtools = "bedtools"
bedgraphtobigwig = "bedGraphToBigWig"

unsorted = sys.argv[1]
genomefile = sys.argv[2]
prefix = sys.argv[3]

#Sort BAM File and create BAM Index file

command1 = samtools+" sort "+unsorted+" "+prefix
cm.getoutput(command1)

#command1_2 = samtools+" view -bq 30 "+prefix+".bam > "+prefix+".filtered.bam"
#cm.getoutput(command1_2)

command2 = samtools+" index "+prefix+".bam "+prefix+".bam.bai"
cm.getoutput(command2)

#Genome-wide bigwig generation
command3 = bedtools+" bamtobed -i "+prefix+".bam > TEMP_BED.bed"
cm.getoutput(command3)

command4 = bedtools+" genomecov -i TEMP_BED.bed -g "+genomefile+" -bga > TEMP_BEDGRAPH.bg"
cm.getoutput(command4)

cm.getoutput("rm TEMP_BED.bed")

command5 = bedgraphtobigwig+" TEMP_BEDGRAPH.bg "+genomefile+" "+prefix+".bw"
cm.getoutput(command5)

cm.getoutput("rm TEMP_BEDGRAPH.bg")


