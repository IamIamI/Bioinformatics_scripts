#!/usr/bin/env Rscript

#A random function found on stackoverflow to remove all currently loaded 
#pacakges... need this to prevent plyr/dplyr conflicts
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics",
                      "package:grDevices","package:utils",
                      "package:datasets","package:methods",
                      "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",
                                                  search()))==1,
                                  TRUE,
                                  FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0) for (package in package.list) detach(package, 
                                                                   character.only=TRUE)
}
detachAllPackages()

#Function to check if package is already installed, if not, installs it. 
# After install it and loads it
install_load <- function(Required_Packages) {
  for(package in Required_Packages){
    if (!package %in% installed.packages()) install.packages(package, 
                                                             character.only = TRUE, 
                                                             dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

install_load(c("getopt","parallel"))

## FUNCTIONs
extract.stats5 <- function(id , path, malt.mode){
    ind <- id
    out <- list()
    for (run in malt.mode){
        ed.dis <- read.table(paste( path,run,'/editDistance/',ind ,'_editDistance.txt', sep="" ),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # keep relevant columns only
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste( path,run,'/damageMismatch/',ind ,'_damageMismatch.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste( path,'default','/readDist/',ind ,'_alignmentDist.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## spec <- 'Yersinia_pestis'
            ed.dis.spec <- ed.dis[ spec , ]
            mp.dam.spec <- mp.dam[ spec , ]
            rd.dis.spec <- rd.dis[ spec , ]
            ## get RatioOfDifferences and N (editdistance0-4 and 0-6)
            ## NOTE: Could be shorter if indeed always only 1 Node presented
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) > 1 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
             mp.dam.spec.max <- max(mp.dam.spec[ rowMax ,"C>T_1"] , mp.dam.spec[ rowMax ,"C>T_2"] , mp.dam.spec[ rowMax ,"G>A_20"], mp.dam.spec[ rowMax ,"G>A_19"])
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.max , read.dis.uniq )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.max , read.dis.uniq )
            }
        }
    }
    out2 <- do.call(rbind,out)
    if(length(malt.mode)==2){
        colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd')
        } else {
            colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd')
        }
    return(out2)
}

local({r <- getOption("repos")
    r["CRAN"] <- "https://mirror.dogado.de/cran" 
    options(repos=r)
})

## CODE

## INFO
## This scripts gathers signatures of species presence at nodes interrogated by MALTextract.
## Evidence is plotted in a heatmap for all samples and solely the samples with species specific evidence.
## For all sample-species pairs with evidence a profile-signature-pdf is plotted.
## Please see profilePDF_explained.pdf for details!
## This script is not designed for 'scan' output and come with no warranty.
## Questions/comments >> key@shh.mpg.de

## USAGE
## Input requires outpath of MALTextract -f <def_anc,default,ancient> run and the taxon-of-interest list (node.list) [e.g. MALTextract taxon list ] .
## Rscript ./postprocessing.v5.r -h

## get options, using the spec as defined by the enclosed list.
## we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    "rmaex.out.fld",  "r" , 1, "character", "MALTextract output folder.",
    "maltex.filter",  "m" , 2, "character", "MALTextract filter mode: <default,def_anc>. This script is not designed for 'scan' output. Default: <def_anc>.",
    "threads",  "t" , 1, "numeric", "Max number of cores used.",
    "help"    ,  "h" , 0, "logical", "Print this help.",
    "node.list"   ,  "n" , 1, "character","List (\\n separated) of nodes to be reported on (aka input species/node list used for MALTextract).",
    "out.name"   ,  "o", 1, "character", "Name of the output file."
), byrow=TRUE, ncol=5);
opt = getopt(spec);


## and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

if ( is.null(opt$out.name) ) { out_file_name <- opt$out.name } else { out_file_name <- "HOPS_excell_summary.csv" }
### ARG parsing and sanity checks
## assign args and modify node.vec (tr ' ' '_')
path <- opt$rmaex.out.fld
if (substr(path,nchar(path),nchar(path)) != "/"){path <- paste(path ,"/",sep="")} # add trailing "/" if missing
if(opt$maltex.filter == 'default') {maltex.mode <- 'default'} else {maltex.mode <- c('default','ancient')}

unq.spec <- unique(gsub(" ","_",scan(file=opt$node,sep="\n",what='character'))) # scan nodes, kill ' ', unique is solely sanity control

## START DATA PROCESSING
all.inds <- colnames(as.matrix(read.table(paste(path,'/default/RunSummary.txt',sep=''),sep="\t",header=T,stringsAsFactors=F,row.names=1,check.names=FALSE,comment.char='')))

### Extract MetaData for all Sample-Species Pairs
out.lists <- mclapply(1:length(all.inds), function(j) extract.stats5( all.inds[j],path,maltex.mode ), mc.cores=opt$threads )
data <- do.call(rbind, out.lists)
data <- data.frame(data,stringsAsFactors=F)
data[, c(4:9,11:16) ] = apply(data[ , c(4:9,11:16)], 2, function(x) as.numeric(as.character(x)))

new_columns <- c('HOPS_step','UniqueToReferences','nonStacked','nonDup','TotalAlignment','destacking?','C>T_1','G>A_-1','mean length (sd)')
data[ , new_columns] <- NA
data[, 'HOPS_step' ] <- 0
data[ data[,'def.dr4'] >= 0.9 & !is.na(data[,'def.dr4']) , 'HOPS_step' ] <- 1## Step1: DiffRatio0-4: > 0.9
data[ data[,'def.mapDam'] > 0 & !is.na(data[,'def.mapDam']) , 'HOPS_step' ] <- 2## Step2: Terminal Damage Present
data[ data[,'anc.dr4'] > 0.8 & !is.na(data[,'anc.dr4']) , 'HOPS_step' ] <- 3## Step3: DiffRatio1-4: > 0.8

folder.names <- paste(path,maltex.mode,'/',sep="")
for (i in 1:length(data[,1])){
  rd <- read.table(paste(folder.names[1],'readDist/',data$id[i],'_alignmentDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
  rd <- rd[grep(data$spec[i], rownames(rd)), ]
  topNode <- rownames(rd)
  if (length(rd[,1])>0){
    ## read mapDamage stuff
    dam <- read.table(paste(folder.names[1],'damageMismatch/',data$id[i],'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    dam <- dam[ grep(data$spec[i],rownames(dam)) , ]
    ## length distribution info
    ld <- read.table(paste(folder.names[1],'readDist/',data$id[i],'_readLengthStat.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ld <- paste( round(ld[ data$spec[i] , 'Mean' ],0),' (',round(ld[ data$spec[i] , 'StandardDev' ],3),')',sep="")
    ## destacking on/off?
    ds <- read.table(paste(folder.names[1],'filterInformation/',data$id[i],'_filterTable.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ds <- ds[data$spec[i], "turnedOn?"]
    ## join output and plot table
    rd2 <- lapply(rd[,-c(1,6)], sum)
#	  data2 <- c(topNode,rd[topNode,'Reference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],rd[topNode,'nonStacked'],ds,round(dam[topNode,'C>T_1'],4),round(dam[topNode,'G>A_20'],4),ld)
#   data2 <- cbind( c('Node','Top Reference','all reads','nonDup','readDis','nonStacked','destacking?','C>T_1','G>A_-1','mean length (sd)'), data)
    data[i,c(18:25)] <- t(rbind(as.matrix(rd2),as.matrix(ds),as.matrix(sum(dam[topNode,'C>T_1'])),as.matrix(sum(dam[topNode,'G>A_20'])),as.matrix(unlist(ld))[1]))
  }
}

if (substr(out_file_name,nchar(out_file_name),nchar(out_file_name)) != ".csv"){out_file_name <- paste(out_file_name ,".csv",sep="")}
write.csv(data, paste(out_file_name,sep=","),row.names = FALSE)
