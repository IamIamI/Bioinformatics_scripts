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
