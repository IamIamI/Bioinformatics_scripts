###
# VCF N corrector
###
# Written on 27-06-2022
# For help contact:
# Lesley Sitter: lesleysitter@hotmail.com
# Alexander Herbig: alexander_herbig@eva.mpg.de
###
# This tool is used to correct VCF files generated by unified genotyper based on 
# references containing N's. UG ignores N's and does not create an entry for these
# in the .vcf file. Tools like SNPevaluator that rely on a complete record of all the 
# positions, will generate frameshifts if this happens. To correct for this we 
# add dummy lines in area's where UG did not generate VCF entries.
#
# Cases that are not handled:
#    - Multi reference VCF files, as SNPevaluation cannot handle these due to the 
#      frame shift
#    - N's at the end of a genome, as the script doesn't know how long the genome is
#      there is metadata in the header but this is not a standard
#
# Script usage should be easy, just invoke Rscript, type the script location and
# add the vcf files you want to correct
# Example: Rscript /path/to/VCF_N_corrector.R file1.vcf file2.vcf
#
# The output uses the orignal name and adds "Ncorrected." to it to identify this
# is the corrected file, and it saves the file in the same location as the original
# vcf.
####

# Everyone must hate R with a burning passions sometimes, while also still loving it
# like a rebelious kid who stole money out of your wallet...
# My current anger is geared towards scientific notation of large numbers...
# Stop that... get some help...
options(scipen=999)
# Also there is one annoying data.frame error that is only conserning the rownames 
# which dont affect the output, therefore i'll supress them. 
oldw <- getOption("warn")
options(warn = -1)
# And if you want to turn warnings back on
# options(warn = oldw)

# Make some space for output
cat("\n")

# Get the commandline arguments, i.e. the files to work on
args = commandArgs(trailingOnly=TRUE)

# Cycle through the files one by one
for(file_in in args){
  message(paste("Start processing: ",file_in))
# file_in <- "C:/Users/lesle/Downloads/SNPevaluation_output/High_genotyping/test.vcf"
# file_in <- "C:/Users/lesle/Downloads/SNPevaluation_output/High_genotyping/AGU007.unifiedgenotyper_test.vcf"
  # Grab the file basename and the directory
  file_name_base <- basename(file_in)
  file_name_path <- dirname(file_in)
  
  # We'll just change the name to reflect that Ns have been correct.
  # Since we dont know if the supplied files is gzipped or what the extention is called
  if (grepl( ".gz", file_in, fixed = TRUE) ){
    vcflist <- read.table(gzfile(file_in))
    if (grepl(".vcf", file_in, fixed = TRUE) ){
      file_name_out <- sub(".vcf.gz", ".Ncorrected.vcf", file_name_base, ignore.case = FALSE)
    }else{
      file_name_out <- sub(".gz", ".Ncorrected.vcf", file_name_base, ignore.case = FALSE)
    }
  }else{
    vcflist <- read.table(file_in)
    if (grepl(".vcf", file_in, fixed = TRUE) ){
      file_name_out <- sub(".vcf", ".Ncorrected.vcf", file_name_base, ignore.case = FALSE)
    }else{
      file_name_out <- paste(file_name_base, ".Ncorrected.vcf", sep="")
    }
  }
  new_file_path <- paste(file_name_path,"/",file_name_out, sep = "")
  

  
  # By default comment lines are ignored when importing files into R, but in a vcf
  # There is a bunch of metadata including the column headers that start with a '#'
  # So we'll load the headers seperately as this is easier and doesnt destroy
  # the column separation for the main table
  # Here we'll also handle a gzipped file again just in case
  if (grepl( ".gz", file_in, fixed = TRUE) ){
    header_metadata <- readLines(gzfile(file_in), 
                                 n=(length(readLines(file_in)) - length(vcflist[,1])))
  }else{
    header_metadata <- readLines(file_in, n=(length(readLines(file_in)) - length(vcflist[,1])))
  }
  
  # Handle the scenario where there are N's at the first positions...
  if (vcflist[1,2]>1){
    missing_data <- data.frame(vcflist[1,1],
                               seq(1,(as.numeric(vcflist[1,2])-1),1),
                               ".","N",".",".",".",".","GT","./.")
    vcflist[,6] <-  as.character(vcflist[,6])
    colnames(missing_data) <- colnames(vcflist)
    vcflist <- rbind(missing_data, vcflist, stringsAsFactors=TRUE)
  }
  
  # Easiest way to identify gaps is by subtracting a position from the position down, and if the value is >1
  # you have a gap. Easiest way is just taking the column, ofsetting it by -1, and then subtracting 
  # the coordinate column with this new -1 column. 
  Coord_offset <- vcflist[-1,2]
  Coord_offset[[length(Coord_offset)+1]] <- Coord_offset[[length(Coord_offset)]]

  # We just remove an additional 1 otherwise every sample has a 1, this way every sample has a 0 and
  # samples with a positive number show the exact number of rows missing
  vcflist$verification <- Coord_offset - (vcflist[,2] + 1)
  
  vcflist_final <- vcflist
  
  # Now the fun part, subset the vcflist to only show ones that have a greater than 0 value 
  for (i in 1:nrow(vcflist[which(vcflist$verification>0),])){
    current_row <- vcflist[which(vcflist$verification>0),][i,]
 
    # Create a new data frame that is exactly the number of missing rows, and contains the proper
    # reference genome name, and site number, as well as dummy empty values for the Ns 
    missing_data <- data.frame(current_row[1],
                               seq((as.numeric(current_row[2])+1),
                                   (as.numeric(current_row[2])+as.numeric(current_row$verification)),
                                   1),
                               ".","N",".",".",".",".","GT","./.","0")
    # Since we are doing this itteratively from bottom to top, and since the row numbers
    # should still match the coordinate numbers, we can just take the whole vcflist from start till the 
    # last row before the missing data, and add the fake data of equal size to the missing data
    colnames(missing_data) <- colnames(vcflist_final)
    vcflist_final <- rbind(vcflist_final[0:as.numeric(current_row[2]),], 
                           missing_data, 
                           tail(vcflist_final,
                                n=(length(vcflist_final[,1])-(as.numeric(current_row[2])))))
  }
  
  
  # Output the data and done
  if (exists("header_metadata") && length(header_metadata)>1){
    write(header_metadata, file=new_file_path, append=FALSE)
    write.table(vcflist_final[,-ncol(vcflist_final)], 
                file=new_file_path, 
                append=TRUE,row.names=FALSE,
                col.names=FALSE,quote=FALSE, sep="\t")
  }else{
    write.table(vcflist_final[,-ncol(vcflist_final)], 
                file=new_file_path, 
                append=FALSE,row.names=FALSE,
                col.names=FALSE,quote=FALSE, sep="\t")
  }
  message(paste("Finished processing: ",file_in))
  cat("\n")
}


