# Python scripts
  
# ClonalFrameML_2_Gff
This script takes a directory in which clonalframeML stored data, and uses
the .newick file and the .importation_status files to generate a gff usefull for feature annotation
and filtering. It also shows which samples fall under which node labels for an easier way to look up data. 
This tool only takes one input, which is the clonalframeML output folder. This folder should contain only one dataset. 
If multiple dataset are stored in the same folder, it's easier to just store them in seperate folders, 
otherwise this script would need more options and validation steps which make it less robust and more convoluted.
  
The script is easily run as followed: 
python ClonalFrameML_2_Gff.py -r </path/to/ClonalFrameML_output/>
  
Output will be stored in the same directory as the input, the gff will be called ClonalFrameML.2.gff
  
# VCF_hetrozygous_positions_barplot
This script is intended to run on a folder full of VCF's (if there is only one VCF in the folder that is fine too), and 
report back the frequency of SNPs called with hetrozygous background noise. The intended purpose is to quickly analyze 
background contamination in bacterial samples, and should off course not be used for eukaryotes or anything like such.
  
This is a very basic python script with no userfriendly error handling or options. 
The script is run as followed 
python VCF_hetrozygous_positions_barplot.py </path/to/VCF_folder/>
  
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

An overview of the expected results of each method
![alt text](https://github.com/IamIamI/Bioinformatics_scripts/edit/master/Python/ssLib_masker.jpg)
