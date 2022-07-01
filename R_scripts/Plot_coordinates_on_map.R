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

# Defining function to draw text boxes
# https://stackoverflow.com/questions/30178954/dynamic-data-point-label-positioning-in-ggmap
draw.rects.modified <- function(d,...){
  if(is.null(d$box.color))d$box.color <- NA
  if(is.null(d$fill))d$fill <- "grey95"
  for(i in 1:nrow(d)){
    with(d[i,],{
      grid.rect(gp = gpar(col = box.color, fill = fill,alpha=0.7),
                vp = viewport(x, y, w, h, "cm", c(hjust, vjust=0.25), angle=rot))
    })
  }
  d
}

# Defining function to determine text box borders
# https://stackoverflow.com/questions/30178954/dynamic-data-point-label-positioning-in-ggmap
enlarge.box.modified <- function(d,...){
  if(!"h"%in%names(d))stop("need to have already calculated height and width.")
  calc.borders(within(d,{
    w <- 0.9*w
    h <- 1.1*h
  }))
}


install_load(c("ggmap","ggpubr"))

my_lon <- c(-61.7745,15.496237,)
my_lat <- c(17.0239,44.078822)
my_years <- c(1998,2020)
my_labels <- c("A","B")


my_df <- data.frame(my_lat,my_lon,my_labels,my_years)

my_map.data <- get_map(location = "world", 
                    maptype = 'roadmap', zoom = 1)
my_map.data <- get_map("world")
#lat_lon <- geocode(c("Hikone-shi Japan", "Aichi Japan", "NA", "Japan", "NA", "NA", "NA", "Japan", "Kochi Japan", "Kochi Japan", "Kyoto Japan", "Kyoto Japan", "NA", "United Kingdom", "United Kingdom", "London United Kingdom", "London United Kingdom", "NA", "NA", "United Kingdom", "NA", "NA", "Exeter United Kingdom", "Aberdeen United Kingdom", "NA", "NA", "USA", "NA", "NA", "NA", "NA", "NA", "Blackpool United Kingdom", "Stockton United Kingdom", "Nottingham United Kingdom", "NA", "NA", "Telford United Kingdom", "NA", "USA", "NA", "United Kingdom", "NA", "London United Kingdom", "London United Kingdom", "NA", "London United Kingdom", "New York USA", "Leiden Netherlands", "London United Kingdom", "NA", "Texas USA", "Memphis USA", "NA", "NA", "NA", "Wales United Kingdom", "London United Kingdom", "United Kingdom", "London United Kingdom", "London United Kingdom", "Wales United Kingdom", "Salisbury United Kingdom", "Alabama USA", "Wales United Kingdom", "United Kingdom", "United Kingdom", "NA", "NA", "France", "Ontario USA", "Gold Coast Australia", "Hong Kong", "Sydney Australia", "Sydney Australia", "NA", "NA", "Wales United Kingdom", "Hong Kong China", "USA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "Canada", "USA", "NA", "NA", "NA", "Brittany France", "Brittany France", "NA", "NA", "NA", "Prague Czechia", "Japan", "Ontario Canada", "NA", "Hong Kong", "NA", "NA", "Georgia USA", "NA", "NA", "NA", "NA", "NA", "NA", "Alencon France", "Italy", "Texas USA", "Hong Kong", "DC USA", "North Carolina USA", "USA", "St.Louis USA", "Brittany France", "Brittany France", "Northern Territory Australia", "Russia", "Russia", "France", "France", "Israel", "Houston USA", "Ontario Canada", "Ontario Canada", "Finland", "Ontario Canada", "Ontario Canada", "Finland", "Finland", "Finland", "Minnesota USA", "Georgia USA", "New Mexico USA", "Connecticut USA", "Minnesota USA", "New York USA", "Maryland USA", "New York USA", "Minnesota USA", "New Mexico USA", "Oregon USA", "Maryland USA", "New Mexico USA", "New Mexico USA", "Finland", "Finland", "Finland", "Finland", "NA", "Christchurch New Zealand", "China", "Manchester United Kingdom", "NA", "Gold Coast Australia", "NA", "NA", "Houston USA", "Houston USA", "Houston USA", "Houston USA", "Houston USA", "Melbourne Australia", "North Carolina USA", "North Carolina USA", "New Zealand", "Australia", "Fiji", "Fiji", "Australia", "Fiji", "Australia", "Fiji", "Kenya", "Brazil", "Fiji", "New Zealand", "Australia", "Fiji", "NA", "India", "Australia", "Australia", "India", "New Zealand", "Fiji", "Brazil", "Australia", "Kenya", "New Zealand", "Fiji", "Australia", "Australia", "Kenya", "New Zealand", "Fiji", "France", "France", "Washington DC USA", "USA", "Brazil", "Brazil", "Australia", "Australia", "Chandigarh India", "Houston USA", "Chandigarh India", "USA", "England United Kingdom", "England United Kingdom", "England United Kingdom", "England United Kingdom", "Australia", "Scotland United Kingdom", "Scotland United Kingdom", "Scotland United Kingdom", "Australia", "Australia"), output = "latlona")
#write.csv(lat_lon, "write_csv.")

boxes <-
  list("top.bumptwice", "calc.boxes",  "enlarge.box.modified", "draw.rects.modified")

# Use GGMAP to plot the samples onto a map
ggmap(get_stamenmap( maptype = c("toner-lite"), bbox = c(left = -180, bottom = -80, right = 179.9999, top = 85), zoom = 3)) +
  geom_point(data = my_df, aes(x = my_lon, y = my_lat, fill = my_years), 
             pch = 21, size = 6) +
  scale_fill_gradient(low="red", high="blue") +
  labs(x = 'Longitude', y = 'Latitude')


#Use geom_density to plot the age distribution onto a graph so that the line smooths out over the number of samples
ggdensity(my_df, x = "my_years", y = "..count..",
          fill = "#0073C2FF", color = "#0073C2FF", rug = TRUE, 
          xlab="Time in years BCE / CE" , ylab = "Frequency")
