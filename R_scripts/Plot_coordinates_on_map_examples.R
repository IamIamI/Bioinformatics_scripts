install_load <- function(Required_Packages) {
  for(package in Required_Packages){
    if (!package %in% installed.packages()) install.packages(package, character.only = TRUE)
    library(package, character.only = TRUE)
  }
}


install_load(c("ggplot2","colorspace","tmaptools","dplyr","maps","tmap","sf","sp"))
# If you get an error that says "there is no package called 'XML'" then run the line below first and rerun the install_load line after
# install.packages("XML", type = "binary")


# Download Spain shapefiles (this probably has all the provinces and cities, keep what you wanna keep and remove what you want to remove)
# https://www.diva-gis.org/gdata
# Download administrative Spain map to get all the regions... And then use the ESP_adm2 shape file which seems to be the province level one
# Lower and higher levels exist, check the .csv's which you want

# Now pick the Spanish .shp file
Spain_file <- file.choose()
#Check if the files actually exists, although this is not really functional, more for my debugging purposes
if (file.exists(Spain_file)){
  print("File accepted")
  Spain <- read_sf(Spain_file)
}

######
# Plotting the empty map
######
tm_shape(Spain) + # Plot the main regions
  tm_layout(bg.color="#FFFFFF", 
            frame="#999999", 
            frame.lwd=3, 
            asp = 1) +
  tm_fill("#E8E8E8") +  # This is the color of the main map
  tm_borders(col = "#999999",  #Province borders around the provinces
             lwd = 1, 
             lty = "solid", 
             alpha = 0.4)


#####
# Remove annoying regions from map
#####
# when looking in the Shapefile or the excell file we can see the province names are mentioned in
# column called "NAME_2". Since these islands are all confined in "NAME_1" as Islas canarias and baleares, we 
# can also remove at that level. 
# We filter out regions we dont care about (weird island or super tiny provinces) by just doing the following
Spain_subsetted <- Spain %>% 
  filter(!NAME_1  %in% c("Islas Canarias", "Islas Baleares"))

tm_shape(Spain_subsetted) + # Plot the main regions
  tm_layout(bg.color="#FFFFFF", 
            frame="#999999", 
            frame.lwd=3, 
            asp = 1) +
  tm_fill("#E8E8E8") +  # This is the color of the main map
  tm_borders(col = "#999999",  #Province borders around the provinces
             lwd = 1, 
             lty = "solid", 
             alpha = 0.4)



#####
# Plotting points on a map using "centroids" i.e. middle points
#####

# This is a HACKY way of copying the shapefile, and changing the coordinates
# The normal centroid approach just puts the coordinate in the middle of each province
Sample_coordinates <- Spain_subsetted %>% 
  filter(NAME_2 %in% c("Granada","Lleida","Almería","Ourense"))
Sample_coordinates_centroids <- st_centroid(Sample_coordinates)
# Now to give weight to the regions... at the plotting stage later on we can
# use this to generate dots or barplots or whatever at the coordinates given by the centroids
# or that we manually changed
Sample_coordinates_centroids$Samples <- c(3,1,2,5)

tm_shape(Spain_subsetted) + # Plot the main regions
  tm_layout(bg.color="#FFFFFF", 
            frame="#999999", 
            frame.lwd=3, 
            asp = 1) +
  tm_fill("#E8E8E8") +  # This is the color of the main map
  tm_borders(col = "#999999",  #Province borders around the provinces
             lwd = 1, 
             lty = "solid", 
             alpha = 0.4) +
tm_shape(Sample_coordinates_centroids) + #These are the positions of the datapoints
  tm_symbols(col = "#b77797", # we can plot all sorts of things, in this case we'll do symbols which is circles by default
             size = "Samples",
             scale = 1.4,
             title.size="Samples per region",
             border.lwd=NA)



#####
# Ofsetting the coordinates of our circles
#####

# To change the coordinates of a given region or in this case a datapoint, we can do the following
# [[1]] signifies it's the first entry, and [1]/[2] signify long and lat or lat, long i forgot which is first
# For example moving the Almería centroid from (-2.3 , 37.2) to (-1.3 , 35)
Sample_coordinates_centroids$geometry[[1]][1] <- -1.3
Sample_coordinates_centroids$geometry[[1]][2] <- 35
# And to change it's value from 3 to 10
Sample_coordinates_centroids$Samples[1] <- 10

tm_shape(Spain_subsetted) + # Plot the main regions
  tm_layout(bg.color="#FFFFFF", 
            frame="#999999", 
            frame.lwd=3, 
            asp = 1) +
  tm_fill("#E8E8E8") +  # This is the color of the main map
  tm_borders(col = "#999999",  #Province borders around the provinces
             lwd = 1, 
             lty = "solid", 
             alpha = 0.4) +
  tm_shape(Sample_coordinates_centroids) + #These are the positions of the datapoints
  tm_symbols(col = "#b77797", # we can plot all sorts of things, in this case we'll do symbols which is circles by default
             size = "Samples",
             scale = 1.4,
             title.size="Samples per region",
             border.lwd=NA)


#####
# Plot a heatmap instead of coordinate based elements
#####

Spain_heatmap <- Spain_subsetted
# i'll just give each province a randomized number between 1 and 100
Spain_heatmap$Heatmap_values <- sample(1:100, 49, replace=TRUE)

tm_shape(Spain_heatmap) + # Plot the main regions
  tm_layout(bg.color="#FFFFFF", 
            frame="#999999", 
            frame.lwd=3, 
            asp = 1) +
  tm_fill(col = "Heatmap_values") +  # This is the color of the main map
  tm_borders(col = "#999999",  #Province borders around the provinces
             lwd = 1, 
             lty = "solid", 
             alpha = 0.4) 