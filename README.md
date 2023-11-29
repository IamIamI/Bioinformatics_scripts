# Bioinformatics_scripts
Scripts made in between projects

# R scripts 
## Merge_files.R

This is a small script that is meant for merging hops summary_table output. 
The use case would be that you might have >100 samples to analyze and to speed up 
the process, you are doing 2 runs of 50 samples each. This results in 2 files.

The script is easy to use, just copy all the RunSummary.txt files from the HOPS 
output, into the same director. And give that directory to the script. 
just invoke it as following
```Rscript Merge_files.R /path/to/your/directory/```


## VCF_N_Corrector.R

This tool is used to correct VCF files generated by unified genotyper based on 
references containing N's. UG ignores N's and does not create an entry for these
in the .vcf file. Tools like SNPevaluator that rely on a complete record of all the 
positions, will generate frameshifts if this happens. To correct for this we 
add dummy lines in area's where UG did not generate VCF entries.

Cases that are not handled:
   - Multi reference VCF files, as SNPevaluation cannot handle these due to the 
     frame shift
    - N's at the end of a genome, as the script doesn't know how long the genome is
     there is metadata in the header but this is not a standard

Script usage should be easy, just invoke Rscript, type the script location and
add the vcf files you want to correct
Example: ```Rscript /path/to/VCF_N_corrector.R file1.vcf file2.vcf fileN.vcf```

The output uses the orignal name and adds "Ncorrected." to it to identify this
is the corrected file, and it saves the file in the same location as the original
vcf.

## Plot_coordinates_on_map.R
A bare bones map for quick coordinate mapping, does not use any google map data etc, making it easy to use

Use at your own perril, is only meant for small data, quick visualization, dirty plotting, and only meant to be used in R Studio with manual editing, no automation provided. 

For nicer looking geo data plots, please refer to 
https://www.lesleysitter.com/2019/08/22/plotting-geo-data/

the code for which is stored on 
https://github.com/IamIamI/pADAP_project/tree/master/geo_plotting_samples 

# Python scripts
## ClonalFrameML_2_Gff
This script takes a directory in which clonalframeML stored data, and uses
the .newick file and the .importation_status files to generate a gff usefull for feature annotation
and filtering. It also shows which samples fall under which node labels for an easier way to look up data. 
This tool only takes one input, which is the clonalframeML output folder. This folder should contain only one dataset. 
If multiple dataset are stored in the same folder, it's easier to just store them in seperate folders, 
otherwise this script would need more options and validation steps which make it less robust and more convoluted.

The script is easily run as followed: 
```python ClonalFrameML_2_Gff.py /path/to/ClonalFrameML_output/```

Output will be stored in the same directory as the input, the gff will be called ClonalFrameML.2.gff

## VCF_hetrozygous_positions_barplot
This script is intended to run on a folder full of VCF's (if there is only one VCF in the folder that is fine too), and 
report back the frequency of SNPs called with hetrozygous background noise. The intended purpose is to quickly analyze 
background contamination in bacterial samples, and should off course not be used for eukaryotes or anything like such.

This is a very basic python script with no userfriendly error handling or options. 
The script is run as followed 
```python VCF_hetrozygous_positions_barplot.py </path/to/VCF_folder/>```

Output will be a file called "VCF_genotyped hetrozygous_loci_frequency.tsv" and will be stored in the same directory as where the script is run from. 

# SSLib_Masker
The script is written for Python 3 and uses the biopython and pysam libraries which can easily be installed by running ```pip install biopython``` and ```pip install pysam``` or using a conda environment.  

This script can mask forward strand 'T' and reverse strand 'A' in Sam/Bam files. The intended goal is to mask ancient "damage" (deaminated cytosines) from singe stranded (SS) libraries and prevent damage from appearing as biological genotypes.  
  
Since the single stranded libraries are not synthetically amplified yet, they are assumed to not have artefactually complemented 'C'>'T'>'A' changes, and instead only have natural deamination artefacts. This means it's possible to mask elements that appear as damage based on the strand that exhibits it.  
  
There are three approaches, the most logical is to mask all the 'T' on the forward strand and the 'A' on the reverse strand as we cannot determine which 'T' are biological and which artefactual.  
The second approach is a reference guide approach where we also supply the reference genome to which the reads were mapped, and we only mask 'T' on the forward strand if the reference sequence has a 'C', and mask the 'A' on the reverse strand if the reference sequence has a 'G'. Although here we assume to know which nucleotides are damage, and will also delete biological C>T, the forward strand 'A' and reverse strand 'T' will still be able to corroborate the genotype if it's biological. Although we do expect that 'C'>'T'/'G'>'A' sites will drop substatially in coverage compared to all other sides, and with low coverage samples might result in a bias during SNP calling.  
  
Lastly we have an edge approach where only a fixed number of nucleotides from the edges are hardmasked. Since deamination apears to mostly happen on the 5' and in a lesser extend on the 3', we could mask all the 'T' on the forwards strand if they are within the first and last couple of nucleotides, and 'A' if they are within the first and last couple of nucleotides of the reverse strand. The user can control this value, and this approach is maybe the most robust. One could use a damage profile as a guide to determine what value to set.  
  
Masking just replaces the nucleotide with an 'N' which does not affect GATK's UnifiedGenotyper when calling genotypes, but might cause problems in other software packages, so use this software at your own discression. 
  
This script has error handling for most thinkable scenario's and should be relatively easy to run as followed:
- Use case: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else  
Example: ```python SSLib_Masker.py --input_file Sample.processed.bam --masking R --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam```  
  
- Use case: When setting the script to Reference guided masking, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position  
Example: ```python SSLib_Masker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam```  
  
- Use case: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.  
Example: ```python SSLib_Masker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam```  
  
- Use case: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0)  
Example: ```python SSLib_Masker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam --mapq_cutoff 37 --len_cutoff 25```  
  
The sofware has an overview of all options which can be called upon by typing 'python SSLib_Masker.py -h' or 'python SSLib_Masker.py --help'  
  
An overview of the options are as followed:  
```
options:
  -m , --masking       Change masking behaviour. 'R' for Reference based Masking. 'H' for HardMasking. 'E' for EdgeMasking. (default: Hardmasking)
  -i , --input_file    Input BAM or SAM file (mandatory)
  -r , --ref_file      Input reference genome file in FASTA format (mandatory for --masking 'R')
  -e , --edge_count    Number of 5' edges to be masked if --masking 'E' is turned on (default: 5)
  -q , --mapq_cutoff   MAPQ cutoff value (default: 0)
  -l , --len_cutoff    Ignore reads below a certain length (default: 0)
  -o , --output_file   Output SAM file with modified reads (default: 'output_modified.sam')
  -h, --help           Show this help message and exit.
```
  
  
# BASH scripts 
## Exp_Desig_Extractor.sh
This script is intended to quickly analyze a BAM file to see which features have coverage and which not. \
This can be very usefull to determine what probes were used for example. The use case was a dataset \
generated by another group, but published only in .BAM format, and in order to compare apples with apples \
we needed to know if our probeset contained the same targets.  
  
To run you require a .BAM file and a .GFF that corresponds with the reference genome that they mapped to   
(if you don't know the reference, you can use 'samtools coverage' for example to see what headers were used, 
and extrapolate from there with a RefSeq/EBI search).\
You also need to have samtools coverage installed / in your PATH to work. It needs to have the samtools coverage\
functionality to work, which was introduced in 1.12 i think?
  
Usage as followed:\
Exp_Desig_Extractor.sh [-g -b -c... -h ] [-g path/to/file ] [-c int]...\
This tool uses Samtools coverage, to determine which gff features were present in any given dataset\
based on BAM alignment.\
&emsp;	-g&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;			Gff-file input [path/to/file]\
&emsp;	-b&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;		Bam file input [path/to/file]\
&emsp;	-c&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;			Coverage threshold, regions below this number will be ignored\
&emsp;	--genesonly&emsp;&emsp;	 Only look at gene features, ignore all other features\
&emsp;	-h | --help&emsp;&emsp;&emsp;		Display help
