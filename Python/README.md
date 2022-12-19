# Python scripts

# ClonalFrameML_2_Gff
This script takes a directory in which clonalframeML stored data, and uses
the .newick file and the .importation_status files to generate a gff usefull for feature annotation
and filtering. It also shows which samples fall under which node labels for an easier way to look up data. 
This tool only takes one input, which is the clonalframeML output folder. This folder should contain only one dataset. 
If multiple dataset are stored in the same folder, it's easier to just store them in seperate folders, 
otherwise this script would need more options and validation steps which make it less robust and more convoluted.

The script is easily run as followed: 
python ClonalFrameML_2_Gff.py /path/to/ClonalFrameML_output/

Output will be stored in the same directory as the input, the gff will be called ClonalFrameML.2.gff
