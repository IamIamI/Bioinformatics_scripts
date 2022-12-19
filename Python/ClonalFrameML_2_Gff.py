#!/usr/bin/python
'''
Program_name: ClonalFrameML_2_Gff

Created on 19-DEC-2022

@author: Lesley Sitter

Description: This script takes a directory in which clonalframeML stored data, and uses
the .newick file and the importation_status files to generate a gff usefull for feature annotation
and filtering. This tool does require the clonalframeML folder to only contain data for one dataset. 
If multiple dataset are stored in the same folder, it's easier to just store them in seperate folders, 
otherwise this script would need more options and validation steps which make it less robust. 
'''

import os
import sys
import glob
import optparse

# Add options here, the names are self explanatory.
def option_parsing(argv):
	usage = ("%prog -r [ClonalFrameML_folder]")
	
	version = ("1.0")
	
	description = ("ClonalFrameML_2_Gff takes the newick file and importation_status files,"
					"and generates a gff containing the 'recombinant positions' with the"
					"corresponding sample to which it applies, plus all names of the"
					"lower level genomes in case positions were identified on a node level")
	
	epilog = ("T.L. Sitter -- lesleysitter@hotmail.com "
				"-- https://nl.linkedin.com/in/lesleysitter")
	
	parser = optparse.OptionParser(usage = usage, description = description, 
									version = "%prog " + version, epilog = epilog)
	
	parser.add_option('-r', '--reads_folder',
		dest = ("folder"),
		metavar = ("FOLDER"),
		help = ("Folder containing the ClonaFrameML output"),
		default = None
		)
	options, remainder = parser.parse_args()
	return(options, remainder)

def option_checker(options):
	# Handle the case if users added a backslash to their directory so we can 
	# consistently deal with this
	if options.folder.endswith("/"):
		pass
	else:
		options.folder = options.folder + "/"
	# First check if the folder the user provided even existst
	if os.path.exists(options.folder):
		# If the folder exists, check if the ClonalFrameML output files we need are there
		file_in_imp_stat = glob.glob(options.folder + "*importation_status.txt")[0]
		try:
			with open(file_in_imp_stat) as file:
				pass
		except:
			print ("\nFile with name %s is not found in folder %s , "
					"please make sure you entered the correct ClonalFrameML output "
					"folder with the option -r \n" % (file_in_imp_stat, options.folder))
			sys.exit(1)
		file_in_newick = glob.glob(options.folder + "*.labelled_tree.newick")[0]
		try:
			with open(file_in_newick) as file:
				pass
		except:
			print ("\nFile with name %s is not found in folder %s , "
					"please make sure you entered the correct ClonalFrameML output "
					"folder with the option -r \n" % (file_in_newick, options.folder))
			sys.exit(1)
	else:
		print ("\nThe folder you specified with option -r '%s' appears to not exists "
				"Please correct and try again\n" % (options.folder))
		sys.exit(1)
	# If none of the check threw an error we should have the two paths to the files we need
	# and we can now return these to variables in the main function
	return(file_in_newick,file_in_imp_stat)

# Stolen from stackoverflow, add keys to dictionary
def add_element(dict, key, value):
	if key not in dict:
		dict[key] = []
	dict[key].append(value)

# Function that process the newick files
def process_newick(file_in_newick):
	node_dict = {}
#	try:
	# Open the output and change the headers
	file_open = open(file_in_newick, 'r').readline()
	# For each nested phylog group (Sample:0.1, Sample:0.1)NODE_01
	# we'll track group position by just checking open and close bracket ()
	count = -1
	# We need to store the sample names temporarily
	sample_label = ""
	# We want to keep track of all the nested positions within our 
	# current node group
	list_pos = []
	# We want to keep track of whether or not we are looking for 
	# opening groups, or closing groups
	node_trigger = False
	label_trigger = False
	# We'll have to go through the newick per character
	for char in file_open:
		# If the character is '(', that means a new phylo group is opening
		if char == "(":
			count += 1
			list_pos.append(count)
			#node_dict[count]=[]
			node_trigger = False
			label_trigger = True
		# Char ':' indicates bootstrap or branchlenght, indicating end of sample names
		elif char == ":":
			# If this name belonged to a node label ')Node_01' we'll change the dictionary
			# kay to this name, so we can easily lookup later which samples are in this node
			if node_trigger == True:
				print (node_dict.keys())
				print(list_pos[-1:][0])
				print(sample_label)
				node_dict[sample_label] = node_dict[list_pos[-1:][0]]
				node_dict.pop(list_pos[-1:][0])
				# We also now have to remove this group from our numerical list
				list_pos = list_pos[:-1]
				label_trigger = False
				node_trigger = False
				sample_label = ""
			# If this name is not a node label, we'll just add the sample name to all the 
			# groups this sample is nested within, higher nodes encompass all lower branched
			# samples
			else:
				for x in list_pos:
					add_element(node_dict,x,sample_label)
				label_trigger = False
				sample_label = ""
		# If there is a ',' that indicates the left side of the pairing is done and a new 
		# sample is coming
		elif char == ",":
			label_trigger = True
			node_trigger = False
		# A closing ')' indicates the end of a group. These are directly followed by a node label
		elif char == ")":
			node_trigger = True
			label_trigger = True
		# If non of the previous characters are present, and the label_trigger is on true, there
		# should be a sample name coming. We will store every character of it in a temp variable
		else:
			if label_trigger == True:
				sample_label = (sample_label + char)
#	except:
#		print ("\nUnable to process ClonalFrameML's newick file. Check if file exists and"
#				"and if it's in it's original format, editing of this file will result in"
#				"incompatibility with this script\n")
#		sys.exit(1) 
	return(node_dict)
	
def process_imp_stat(file_in_imp_stat, options, node_dict):
	file_open = open(file_in_imp_stat, 'r').readlines()
	file_out = open(options.folder + "ClonalFrameML.2.gff", 'w')
	file_out.write("##gff-version 3\n")
	file_out.write("##sequence-region SEQUENCE 1 Unknown\n")
	for lines in file_open:
		if "Node\tBeg\tEnd" in lines:
			pass
		else:
			line = lines.split('\t')
			if "NODE" in lines: 
				print (lines)
				print (node_dict[line[0]])
				gff_line = ("SEQUENCE\tClonalFrameML\tCDS\t" + line[1] + "\t" + line[2].strip() + "0.000\t.\t0\tID=\"" + line[0] + ":" + ','.join(node_dict[line[0]]) + "\"\n")
				file_out.write(gff_line)
			else:
				gff_line = ("SEQUENCE\tClonalFrameML\tCDS\t" + line[1] + "\t" + line[2].strip() + "0.000\t.\t0\tID=\"" + line[0] + "\"\n")
				file_out.write(gff_line)
	
def main():
#	# Check which Python version is being run
#	cur_version = sys.version_info
#	# If version <3 is ran, options might need to change
#	if cur_version < (2,9):
#		options.pyversion = 2
#	else:
#		options.pyversion = 3
		
		
	# Get the options
	options, remainder = option_parsing(sys.argv[1:])
	# Check correctness of options, and if files exists etc
	# Return the newick file and important_status file paths
	file_in_newick, file_in_imp_stat = option_checker(options)
	
	# Submit the newick file, extract node labels, find which 
	# samples fall under that label and return them as a dictionary
	node_dict = process_newick(file_in_newick)
	
	# Process important statistic files, also add the node_dict 
	# so that we can replace NODE labels with the sample names
	process_imp_stat(file_in_imp_stat, options, node_dict)
	
	print ("\nOutput files should now be available in %s.\n"
			"The output file is called 'clonalframeml.2.gff'" % (options.folder))
	sys.exit(1)

if __name__ == "__main__":
	main()
