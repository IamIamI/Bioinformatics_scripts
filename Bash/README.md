#BASH scripts 
#Exp_Desig_Extractor.sh
This script is intended to quickly analyze a BAM file to see which features have coverage and which not.
This can be very usefull to determine what probes were used for example. The use case was a dataset 
generated by another group, but published only in .BAM format, and in order to compare apples with apples
we needed to know if our probeset contained the same targets.

To run you require a .BAM file and a .GFF that corresponds with the reference genome that they mapped to 
(if you don't know the reference, you can use 'samtools coverage' for example to see what headers were used, and extrapolate from there
with a RefSeq/EBI search)
You also need to have samtools coverage installed / in your PATH to work. It needs to have the samtools coverage
functionality to work, which was introduced in 1.12 i think? 

Usage as followed:
Usage: 
Exp_Desig_Extractor.sh [-g -b -c... -h ] [-g path/to/file ] [-c int]...
This tool uses Samtools coverage, to determine which gff features were present in any given dataset
based on BAM alignment. 
	-g			Gff-file input [path/to/file]
	-b			Bam file input [path/to/file]
	-c			Coverage threshold, regions below this number will be ignored
	--genesonly		Only look at gene features, ignore all other features
	-h | --help		Display help