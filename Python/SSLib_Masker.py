#!/usr/bin/python

# Script name: SSLib_Masker.py
# Version: v1.0
# Author: Thomas Lesley Sitter
# Contact: lesleysitter@hotmail.compile

# Description: Used to mask thymine and adenine sites in ancient SAM/BAM files created from single stranded libraries, to prevent damage from being incorporated into genotyping.

# Usability: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else
# Example: python SSLib_Masker.py --input_file Sample.processed.bam --masking R --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam
# Usability: When setting the script to Reference guided masking, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position
# Example: pythonSSLib_Masker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam
# Usability: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.
# Example: python SSLib_Masker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam
# Optional: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0)
# Example: python SSLib_Masker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam --mapq_cutoff 37 --len_cutoff 25

# Comment: Use at your own discression, created for personal use
#

import os
import sys
import math

# Check if the script is running with Python 3, this first one seems to not work?? 
if sys.version_info.major < 3:
	print("\nError: This script requires Python 3. Please run it with a Python 3 interpreter.\n\n")
	sys.exit(1)
# Check if the script is running with Python 3
if not sys.version_info[:1] == (3):
	print("\nError: This script requires Python 3. Please run it with a Python 3 interpreter.\n\n")
	sys.exit(1)

import argparse

# Check if the pysam library is installed
try:
	import pysam
except ImportError:
	print("\nError: This script requires pysam to be installed.\nTry \'pip install pysam\' or \'pip3 install pysam\' (depending on your setup) to install it.\n\n")
	sys.exit(1)

# Check if the pysam library is installed
try:
	from Bio import SeqIO
except ImportError:
	print("\nError: This script requires biopython to be installed.\nTry \'pip install biopython\' or \'pip3 install biopython\' (depending on your setup) to install it.\n\n")
	sys.exit(1)

# Determine if a bam or sam is provided
def sam_bam_picker(input_file,ref_file,output_file,mapq_cutoff,len_cutoff,masking,edge_count):

	type_testing = input_file.lower()
	# Check if the read is mapped and MAPQ is above the cutoff
	if '.sam' in type_testing :
		with pysam.AlignmentFile(input_file, 'r') as sam, pysam.AlignmentFile(output_file, 'w', header=sam.header) as output_sam:
			process_sam_bam(sam, output_sam, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count)

	elif '.bam' in type_testing :
		with pysam.AlignmentFile(input_file, 'rb') as bam, pysam.AlignmentFile(output_file, 'wb', header=bam.header) as output_bam:
			process_sam_bam(bam, output_bam, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count)

	else:
		print(f"\nError: It could not be determined if '{input_file}' is a SAM or BAM formatted file. Make sure the file has a .bam or .sam extention.\n\n")
		sys.exit(1)

def process_sam_bam(in_file, out_file, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count):
	if masking == "R":
		# Load the reference genome
		reference_dict = SeqIO.index(ref_file, "fasta")

	# Go through the reads in the FASTA
	for read in in_file:
		if not read.is_unmapped and read.mapping_quality >= mapq_cutoff and read.query_length >= len_cutoff:
			read_sequence = read.query_sequence # Get read sequence
			modified_qualities = read.query_qualities # Obtain the original quality values

			if masking == "R":
				# Load the reference sequence
				reference_seq = reference_dict[read.reference_name][read.reference_start:read.reference_start+len(read_sequence)]

				# Forward read
				if not read.is_reverse:
					modified_sequence = ''.join(['N' if ref == 'C' and read == 'T' else read for ref, read in zip(reference_seq, read_sequence)])
				# Reverse read
				else:
					modified_sequence = ''.join(['N' if ref == 'G' and read == 'A' else read for ref, read in zip(reference_seq, read_sequence)])

			elif masking == "H":
				# Forward read
				if not read.is_reverse:
					modified_sequence = ''.join(['N' if read == 'T' else read for read in read_sequence])
				# Reverse read
				else:
					modified_sequence = ''.join(['N' if read == 'A' else read for read in read_sequence])

			elif masking == "E":
				if edge_count >= math.floor(len(read_sequence)/2):
					mask_numb = math.floor(len(read_sequence)/2)-1
				else:
					mask_numb = edge_count
				# Forward read
				if not read.is_reverse:
					modified_seq_L = ''.join(['N' if read == 'T' else read for read in read_sequence[:mask_numb]])
					modified_seq_R = ''.join(['N' if read == 'T' else read for read in read_sequence[len(read_sequence)-mask_numb:]])
					modified_sequence = modified_seq_L + read_sequence[mask_numb:len(read_sequence)-mask_numb] + modified_seq_R
				# Reverse read
				else:
					modified_seq_L = ''.join(['N' if read == 'A' else read for read in read_sequence[:mask_numb]])
					modified_seq_R = ''.join(['N' if read == 'A' else read for read in read_sequence[len(read_sequence)-mask_numb:]])
					modified_sequence = modified_seq_L + read_sequence[mask_numb:len(read_sequence)-mask_numb] + modified_seq_R

			else:
				print(f"\nError: The masking setting '{masking}' is not recognized. Please use \'S\' for SoftMasking, \'H\' for HardMasking, and \'E\' for EdgeMasking.\n\n")
				sys.exit(1)

			read.query_sequence = modified_sequence  # Update the sequence field
			read.query_qualities = modified_qualities # Preserve the original quality values
			
			out_file.write(read)

def main():
	parser = argparse.ArgumentParser(description='Mask a SAM/BAM file for deaminated bases based on reference genome. The script can softmask, hardmask and edgemask', add_help=False)
	parser.add_argument('-m', '--masking', default="H", metavar='', help='Change masking behaviour.\n\'R\' for Reference based Masking.\n\'H\' for HardMasking.\n\'E\' for EdgeMasking. (default: Hardmasking)')
	parser.add_argument('-i', '--input_file', metavar='', help='Input BAM or SAM file (mandatory)')
	parser.add_argument('-r', '--ref_file', default="NA", metavar='', help='Input reference genome file in FASTA format (mandatory for --masking \'R\')')
	parser.add_argument('-e', '--edge_count', type=int, default=5, metavar='', help='Number of 5\' edges to be masked if --masking \'E\' is turned on (default: 5)')
	parser.add_argument('-q', '--mapq_cutoff', type=int, default=0, metavar='', help='MAPQ cutoff value (default: 0)')
	parser.add_argument('-l', '--len_cutoff', type=int, default=0, metavar='', help='Ignore reads below a certain length (default: 0)')
	parser.add_argument('-o', '--output_file', metavar='', default='output_modified.sam', help='Output SAM file with modified reads (default: \'output_modified.sam\')')
	parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
	args = parser.parse_args()

	# Check if the input BAM file exists
	if not os.path.exists(args.input_file):
		print(f"\nError: Input SAM/BAM file '{args.input_file}' not found.\n\n")
		return
		# Check if the input BAM file exists

	if args.masking == "R":
		if not args.ref_file =="NA":
			if not os.path.exists(args.ref_file):
				print(f"\nError: Input reference FASTA file '{args.ref_file}' not found.\n\n")
				return
		else:
			print(f"\nError: When using the \'R\' Reference guided masking, the user has to supply a reference genome with --ref_file.\nNo Reference was supplied.\nPlease supply the reference genome against which the reads were mapped,\nor choose the \"H\" or \"E\" masking option instead by setting these using --masking.\n\n")
			return
	
	sam_bam_picker(args.input_file, args.ref_file, args.output_file, args.mapq_cutoff, args.len_cutoff, args.masking, args.edge_count)

if __name__ == "__main__":
	main()
