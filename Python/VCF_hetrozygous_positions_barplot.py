import os
import sys
import glob
import re

# Arg 1 is a directory containing VCF files

vcf_folder = sys.argv[1]

if not vcf_folder.endswith("/"):
    vcf_folder = (vcf_folder + "/")
    
file_out = open("VCF_genotyped hetrozygous_loci_frequency.tsv", 'w')
file_out.write("Sample\t0.0\t0.1\0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8\t0.9\t1.0\n")

# Now for the VCF files we want to analyze
for file_path in glob.glob(vcf_folder + '*/*.unifiedgenotyper.vcf'):
    basename = os.path.basename(file_path)
    basename = re.sub(r'.unifiedgenotyper.vcf$', '', basename)
    
    print("Processing sample: %s" % (basename))
    # We'll reformat the filename/location etc a bit to make sure we don't overwrite the original vcf
    file_open = open(file_path, 'r').readlines()
    coordinate_list=[]
    freq_dict = {0.0:0, 0.1:0, 0.2:0, 0.3:0, 0.4:0, 0.5:0, 0.6:0, 0.7:0, 0.8:0, 0.9:0, 1.0:0}

    
    for lines in file_open:
        if not lines.startswith('#'):
            line = lines.split('\t')
            if "GT:AD" in line[8]:
                loci_freq = line[9].split(":")[1]
                loci_freq = loci_freq.split(",")
                loci_freq = [int(x) for x in loci_freq]
                loci_freq = round((loci_freq[-1]/sum(loci_freq)),1)
                freq_dict[loci_freq] += 1
    file_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (basename,freq_dict[0.0],freq_dict[0.1],freq_dict[0.2],freq_dict[0.3],freq_dict[0.4],freq_dict[0.5],freq_dict[0.6],freq_dict[0.7],freq_dict[0.8],freq_dict[0.9],freq_dict[1.0]))
