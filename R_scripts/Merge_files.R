#A function that can detach packages, use this to resolve conflicts for example between plyr and dplyr
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

#Function to check if package is already installed, if not, installs it. After install it and loads it
install_load <- function(Required_Packages) {
  for(package in Required_Packages){
    if (!package %in% installed.packages()) install.packages(package, character.only = TRUE)
    library(package, character.only = TRUE)
  }
}

#install and load plyr
install_load("plyr")

#Pick the directory
directory <- choose.dir()

#Make a list of all the files in the directory
LF <- list.files(path=directory,pattern=".txt",full.names=T,recursive=T)

#Just store the first LF entry in a data_frame
sdata1 <- read.table(LF[1], sep="\t", header=TRUE)

#Remove the entry from the list so we don't add the first file twice
LF2 <- LF[-1]

#We make a function that just reads in a file and merges it with the previously
#create sdata1
merge_function <- function(x, out){
  x <- read.table(x, sep="\t", header=TRUE)
  y <- merge(sdata1, x, by = "Node", all = TRUE)
  return(y)
}

#We recursively run through the files in LF2 and preform the previously written merge function
for(i in LF2){
  sdata1 <- merge_function(i)
}

#For my particular case i want the data transposed since that will create less colums (easier on the eyes)
sdata1 <- t(sdata1)

#Finally, output the data
write.csv(sdata1, file="Something.csv")
