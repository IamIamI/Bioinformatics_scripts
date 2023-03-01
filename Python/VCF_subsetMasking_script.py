import os
import sys
import glob
import re

# Arg 1 is the Gff you want to use to mask the VCF files
# Arg 2 is a directory containing VCF files
file_open = open(sys.argv[1], 'r').readlines()

def checkIfinLookup(mylist):
    for x in mylist:
        if "taxa=" in x:
            return x

filter_dictionary = {}

# First process the gff files line by line
for lines in file_open:
    if not lines.startswith('#'):
        lines = lines.split('\t')
        # Extract start and end coordiantes of the feature
        start_coord = lines[3]
        end_coord = lines[4]
        # The for Gubbins the metadata is more bulky, we only care about the taxanomic names though
        # so split the data by ';', look for the taxa part, take that and remove all the " quote marks
        if "GUBBINS" in lines:
            meta_data = lines[8].split(';')
            meta_data = checkIfinLookup(meta_data)
            meta_data = re.sub(r'^taxa=\"', '', meta_data)
            meta_data = re.sub(r'\"$', '', meta_data)
            # Finally, samples are separated by uneven amounts of spaces... annoying... 
            # so split by spaces, then remove the entrys that are empty (i.e. if you have 3 spaces,
            # and split by space, you'll get one or two empty entries in your list... remove those
            meta_data = meta_data.split(' ')
            sample_names=([sample for sample in meta_data if len(sample)>1])
        # If you ran ClonalFrameML2_2_Gff.py the formatting is simpler, it only contains a ID=
        # followed by the sample names seperated by comma's. If it concers a NODE with multiple 
        # samples, its written as NODE_##:Sample1,Sample2 ... so we can remove this NODE part by
        # first splitting by : and taking only the back portion... then seperating the rest by comma
        elif "ClonalFrameML" in lines:
            meta_data = re.sub(r'^ID=\"', '', lines[8])
            meta_data = meta_data.strip('"\n')
            if ":" in meta_data:
                meta_data = meta_data.split(':')[1]
                meta_data = meta_data.split(',')
        else:
            continue
        # Lastly now that we have all the formatting the same way (i.e. each line is just a list 
        # of samples, a start and stop coordiante and nothing more. We can now add this to a dictionary
        for sample in meta_data:
            if sample in filter_dictionary :
                filter_dictionary[sample].append("%s-%s" % (start_coord,end_coord))
            elif len(sample) > 1 and not sample.isspace():
                filter_dictionary[sample]=[("%s-%s" % (start_coord,end_coord))]

directory = sys.argv[2]
if not directory.endswith("/"):
    directory = (directory + "/")
    
    
[print(keys) for keys in filter_dictionary]

# Now for the VCF files we want to analyze
for file_path in glob.glob(directory + '*/*.vcf'):
    basename = os.path.basename(file_path)
    basename = re.sub(r'.unifiedgenotyper.vcf$', '', basename)
    
    # if we don't have any filtering parameters for a particular sample, there is no need to even open it
    if basename in filter_dictionary:
        print("Processing sample: %s" % (basename))
        # We'll reformat the filename/location etc a bit to make sure we don't overwrite the original vcf
        copy_file = (file_path.strip(".vcf")+".TwoX_subsetMasked_v2.vcf")
        file_open = open(file_path, 'r').readlines()
        file_out = open(copy_file, 'w')
        mask_coordinates = filter_dictionary[basename]
        for lines in file_open:
            if lines.startswith('#'):
                file_out.write(lines)
            else:
                line = lines.split('\t')
                for entry in mask_coordinates:
                    entry = entry.split('-')
                    start_coord = int(entry[0])
                    end_coord = int(entry[1])
                    if int(line[1]) >= start_coord and int(line[1]) <= end_coord:
                        lines = ("%s\t%s\t%s\t%s\t.\t.\t.\t.\tGT\t./.\n" % (line[0],line[1],line[2],line[3]))
                file_out.write(lines)
                

    else:
        continue