# A random function found on stackoverflow to remove all currently loaded pacakges... need this to prevent plyr/dplyr conflicts
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

#Install and load the following packages
install_load(c("plyr"))
remove.packages("yaml")
install_load(c("yaml"))
#Load a folder of files generated with 
files_list <- list.files(path=choose.dir(), pattern="*.tsv", full.names=TRUE, recursive=FALSE)
read_file <- read.csv(file.choose(), sep="\t", na.strings=c("","NA"), header = FALSE)

coverage_distro <- data.frame(Name=character(),
                              Coverage =character(), 
                              Frequency=character(),
                              stringsAsFactors=TRUE) 

average_coverage_data <- data.frame(Name=character(),
                                    Ref_nuc_span=character(), 
                                    Samp_nucs_mapped=character(),
                                    Average_coverage=character(),
                                    stringsAsFactors=TRUE) 

for (file in files_list) {
  
  table_input <- as.data.frame(read.csv(file, sep="\t", na.strings=c("","NA"), header=FALSE))
  #extract sample numbers
  colnames(table_input)[4] <- "Coverage"
  sample_name = sub("_coverage.tsv","",sub(".*/","",file))
  #column bind, sample name, the coverage value, and the length of the block
  table_input$Name <- sample_name
  table_input$Ref_nuc_span <- c(table_input[,3]-table_input[,2])
  table_input <- table_input[,-c(1:3)]
  #add column containing multiplecation of coverage times nucleotides in block to get total coverage over that area
  table_input$Samp_nucs_mapped <- c(table_input$Coverage*table_input$Ref_nuc_span)
  
  #Make a new table containing the name, total ref nucleotids covered, total sample nucleotides aligned and average coverage
  coverage_data <- data.frame(Name = sample_name,
                              Ref_nuc_span = sum(table_input[,"Ref_nuc_span"]),
                              Samp_nucs_mapped = sum(table_input[,"Samp_nucs_mapped"]),
                              Average_coverage = sum(table_input[,"Samp_nucs_mapped"])/sum(table_input[,"Ref_nuc_span"]))
  
  #join the coverage list to the average_coverage_data table
  average_coverage_data <- rbind(average_coverage_data,coverage_data)
  
  #Join the coverage distribution to an coverage_distro table
  coverage_distro <- rbind(coverage_distro,cbind(sample_name,ddply(table_input, .(Coverage), function(x) sum(x[,3]))))                          
}


detachAllPackages()
install_load(c("ggplot2", 
               "dplyr",
               "cowplot",
               "viridis"))

# Add a column of Reads that we obtained seperately to the average_coverage_data file
colnames(read_file) <- c("Name","Reads")
read_file$Reads <- c(read_file$Reads * 2)
average_coverage_data2 <- merge(average_coverage_data, read_file, by.x="Name", by.y="Name", sort = TRUE)

# Add the mean total ref nucleotids covered, used to make a line in our barchart
average_coverage_data2 <- average_coverage_data2 %>% mutate(Mean_ref_nuc_span = mean(Ref_nuc_span))
average_coverage_data2 <- average_coverage_data2 %>% mutate(Mean_tot_nucs_aligned = mean(Samp_nucs_mapped))

# This is used for our triple bar plot where we compare the total sample reads mapped to the total ref 
# nucleotides covered and total unmapped reads
average_coverage_tmp1 <- average_coverage_data2 %>% mutate(Mean_Samp_nucs_mapped = mean(Samp_nucs_mapped))
average_coverage_tmp1$label <- "Mapped"
average_coverage_tmp2 <- average_coverage_tmp1 %>% mutate(Samp_nucs_mapped = Reads*100)
average_coverage_tmp2$label <- "Est. Total"
average_coverage_tmp3 <- average_coverage_tmp1 %>% mutate(Samp_nucs_mapped = c((Reads*100) - Samp_nucs_mapped ))
average_coverage_tmp3$label <- "Unmapped"
average_coverage_tmp <- rbind(average_coverage_tmp1,average_coverage_tmp2,average_coverage_tmp3)

# These are the average coverage per sample as well as the 10x and 30x indicaters and the mean coverage of all samples
# These values will be used to generate indicater lines in the boxplot
average_coverage_data2 <- average_coverage_data2%>%  mutate(mean_cov = mean(Average_coverage))
average_coverage_data2$tenx = 10
average_coverage_data2$thirtyx = 30

# Lastly, the mean reads, used to generate the indicater line later on
average_coverage_data2 = average_coverage_data2 %>% mutate(mean_rd = mean(Reads))

colnames(coverage_distro) <- c("Name","Coverage","Frequency")
# The first plot generates the coverage distribution per sample, so these will be 100 little plots in one
pdf("Coverage_distributions.pdf", width=18,height=14)
p1 <-   ggplot(coverage_distro, aes(x=Coverage, y=Frequency, group=Name, fill=Name)) +
        geom_area() +
#        scale_fill_viridis(discrete = TRUE) +
        scale_x_continuous(limit=c(0,30)) +
        theme(legend.position="none") +
        theme(legend.position="none",
              panel.spacing = unit(0.1, "lines"),
              strip.text.x = element_text(size = 8)) +
        facet_wrap(~Name, scale="free_y")
p1
dev.off()


pdf("Sample_read_quality.pdf", width=18,height=12)
# The following plots will generate our barcharts and median lines, as well as additional information we added
plot1 <-  average_coverage_data2 %>%
          select(Name, Ref_nuc_span  ) %>%
          na.omit() %>%
          ggplot() +
          geom_bar(aes(x = Name, y = Ref_nuc_span  ), stat = "identity", alpha = 0.75, fill = "palegreen3") +
          ylab("Chrom Nuc covered") +
          geom_errorbar(data=average_coverage_data2, aes(Name, ymax = Mean_ref_nuc_span, ymin = Mean_ref_nuc_span),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          theme_minimal() +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())

plot2 <-  average_coverage_tmp %>%
          select(Name, Samp_nucs_mapped,label) %>%
          na.omit() %>%
          ggplot() +
          geom_bar(aes(x = Name, y = Samp_nucs_mapped , fill=label), stat = "identity", position=position_dodge()) +
          ylab("Read Nucs mapped") +
          geom_errorbar(data=average_coverage_tmp, aes(Name, ymax = Mean_Samp_nucs_mapped, ymin = Mean_Samp_nucs_mapped),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          scale_fill_brewer(palette="Accent")+
          theme_minimal() +
          guides(shape = guide_legend(override.aes = list(size = 100)),
                 color = guide_legend(override.aes = list(size = 100))) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
        #        legend.position = "none")
                legend.title = element_blank(), 
                legend.text=element_text(size=4),
                legend.position = c(0.95, 0.8))

cols <- c("Mean"="black","10X"="maroon","30X"="seagreen1")
plot3 <-  average_coverage_data2 %>%
          select(Name, Average_coverage) %>%
          na.omit() %>%
          ggplot() +
          geom_bar(aes(x = Name, y = Average_coverage), stat = "identity", alpha = 0.75, fill = "slategray3") +
          ylab("Avg coverage") +
          geom_errorbar(data=average_coverage_data2, aes(Name, ymax = mean_cov, ymin = mean_cov, colour="Mean"),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          geom_errorbar(data=average_coverage_data2, aes(Name, ymax = tenx, ymin = tenx, colour="10X"),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          geom_errorbar(data=average_coverage_data2, aes(Name, ymax = thirtyx, ymin = thirtyx, colour="30X"),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          theme_minimal() +
          scale_colour_manual(name="Error Bars",values=cols) + 
          scale_fill_manual(name="Bar",values=cols) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                legend.title = element_blank(), 
                legend.text=element_text(size=4),
                legend.position = c(0.90, 0.8))

plot4 <-  average_coverage_data2 %>%
          select(Name, Reads) %>%
          na.omit() %>%
          ggplot() +
          geom_bar(aes(x = Name, y = Reads), stat = "identity", alpha = 0.75, fill = "palevioletred2") +
          ylab("# Reads") +
          geom_errorbar(data=average_coverage_data2, aes(Name, ymax = mean_rd, ymin = mean_rd),
                        size=0.5, linetype = "longdash", inherit.aes = F, width = 1) +
          theme_minimal() +
          theme(axis.title.x = element_blank())

# Finally we can plot all the 4 plots on top of each other and bask in the glory of R visualizations
plot_grid(plot1, plot2, plot3, plot4, align = "v", ncol = 1, rel_heights = c(0.24, 0.24, 0.24, 0.28))
dev.off()