### Geography 518: Advanced Spatial Analysis ###
#FINAL PROJECT

install.packages("sf")
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("spatstat")
install.packages("sp")
install.packages("tmap")
install.packages("gstat")
install.packages("BAMMtools")
install.packages("shinyjs")
install.packages("spgwr")

library(sf)
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(tmap)
library(gstat)
library(BAMMtools)
library(shinyjs)
library(spgwr)

dir <- "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master"
setwd(dir)


#-----------------------------------------------------------------------------

#------------------------------------- THIS SECTION UPLOADS & CLEANS CENSUS DATA

#Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data

pm25 <- pm25[,1:2] #Select only columns 1 and 2
colnames(pm25) <- c("POSTALCODE", "PM25") #Change the column names 
pm25 <- na.omit(pm25)

#Reading in postal code shapefile
postalcodes <- readOGR("./BC_PostalCodes","BC_Postal_Codes") #Read in related postal code data

#Reading in dissemination tract and income data
income <- read.csv("Income.csv") #Read in census income data  
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns

census.tracts <- readOGR("./BC_DA","BC_DA") #Read in dissemination tract shapefile

income.tracts <- merge(census.tracts,income, by = "DAUID") #Merge income and dissemination area data
nrow(income.tracts) #Determine the number of columns in the dataframe
income.tracts <- income.tracts[!is.na(income.tracts$Income),]

crs(income.tracts)

# Create choropleth map of income
# med.income <- income.tracts$Income
# shades <- auto.shading(med.income, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, med.income, shades) #map the data with associated colours
# choro.legend(-123.5764,49.18482, shades, cex = 0.8) #add a legend (you might need to change the location)
# title("Median Income in Vancouver, British Columbia")


#------------------------------------------------------------------ STUDY SITE MAP

STUDYSITE <- tm_shape(income.tracts) + 
  tm_polygons( 
            title="Greater Vancouver Regional District",
            style = "fisher", #suggested by R for larger datasets
            palette = "Oranges", n = 6, border.col = "black", lwd = 0.3, alpha = 0.0)+
  tm_scale_bar(width = 0.25, position = c("left","bottom"))
tmap_mode("view")
STUDYSITE
#png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/STUDYSITEMAP.png", width = 1000, height=800)
#dev.off()

#----------------------------------

#-----------------------------------Alternative Mapping Technique - INCOME CHLOROPLETH MAP

map_MedIncome <- tm_shape(income.tracts) + 
  tm_polygons(col = "Income", 
              title = "Median Income\nGreater Vancouver Regional District, British Columbia",
              style = "fisher", #suggested by R for larger datasets
              palette = "viridis", n = 6, alpha =0.8, border.alpha = 0)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))

#tmap_mode("view")
tmap_mode("plot")
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/IncomeChoroplethMap.png", width = 1000, height=800)
map_MedIncome #Local Indicator Spatial Autocorrelaetion
dev.off()


#-------------------------------------------------------------

#Select postal codes that fall within dissemination tracts
postalcodes <- intersect(postalcodes,income.tracts)
plot(postalcodes) #See what the data looks like spatially
head(postalcodes) #See what the data looks like in tabular form

#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")

#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the max.
#essentially choosing the max PM2.5 value from each Dissemination Area to represent the DA.   
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")

#Create a subsample of the datapoints provided in the PM2.5 dataset using the sample 'n' provided on CourseSpaces
set.seed(160)
sampleSize=160
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]
minPM25 <- min(spSample$PM25AGG, na.rm=T)
maxPM25 <- max(spSample$PM25AGG, na.rm=T)
#make min & max

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(spSample, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)


#-------------------------------------------------------------------- Plot the sample sites of the map
STUDYSITE_WITH_PTS <- tm_shape(income.tracts) + 
  tm_polygons(title="PM2.5 Sample Sites",
              palette = "Oranges", n = 6, border.col = "black", lwd = 0.3, alpha = 0.0) +
  tm_shape(spSample) + tm_dots(size=0.009) +
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))

tmap_mode("view")
STUDYSITE_WITH_PTS
tmap_mode("plot")
#------------------------------------------------------------------------------- Global & Local Moran's I -----------

#Weighted Matrix, spatially based on Census Tract information

vancouver.nb <- poly2nb(income.pm25) #using study area (defined by income) tracts, queens weighting -> defining a polygon's neighbourhood, up to eight surrounding points 
vancouver.net <- nb2lines(vancouver.nb,coords=coordinates(income.pm25))

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/QueensNeighbourhood.png", width = 800, height=600)
tm_shape(income.pm25) + 
  tm_borders(col='lightgrey') + 
  tm_shape(vancouver.net) + 
  tm_lines(col='blue') #creating a visualization of the Queens case weighted mean
dev.off()


#creating the weight matrix (i.e a bunch of 1's and 0's, defining if they're neighbours or not)
vancouver.lw <- nb2listw(vancouver.nb, zero.policy = TRUE, style = "W")
print.listw(vancouver.lw, zero.policy = TRUE)



#---------------------------------------------------------------- Global Moran's I Test (i.e. Global Spatial Autocorrelation)

mi <- moran.test(income.pm25$Income, vancouver.lw, zero.policy = TRUE) #based off vancouver.lw weighted matrix
mi


moran.range <- function(lw) { 
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vancouver.lw) #we are going to need these values to calculate the Z-score


mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- ((mI-eI)/(sqrt(var)))
z

#------------------------------------------------------------------ Local Moran's I Test (i.e. Local Spatial Autocorrelation)

lisa.test <- localmoran(income.pm25$Income, vancouver.lw) #based off Vancouver.lw weighted matrix

#creating columns for expected values 
income.pm25$Ii <- lisa.test[,1] #this describes the expected values we are going to get: I, E,Var, Z,P
income.pm25$E.Ii<- lisa.test[,2]
income.pm25$Var.Ii<- lisa.test[,3]
income.pm25$Z.Ii<- lisa.test[,4]
income.pm25$P<- lisa.test[,5]

#mapping the local I 

map_LISA <- tm_shape(income.pm25) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I value of Median Income,\nFisher Classification", 
              style = "fisher", #suggested by R for larger datasets
              palette = "-RdBu", n = 5, lwd = 0.5,
              midpoint = NA) +
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))

tmap_mode("plot")
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/LISA_MoransIValue_IncomeMap.png", width = 800, height=600)
map_LISA #Local Indicator Spatial Autocorrelaetion
dev.off()

#plot based on p - values : which indicate clustering (low p value) - define breaks 
map_LISA_p <- tm_shape(income.pm25) + 
  tm_polygons(col = "P", 
              title = "Local Moran's I p-value of Median Income,\nFisher Classification", 
              style = "fisher", #suggested by R for larger datasets
              palette = "-RdBu", n = 5, lwd = 0.5) +
  tm_scale_bar(width = 0.25, position = c("left","bottom"))

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/LISA_PValue_IncomeMap.png", width = 800, height=600)
map_LISA_p #Local Indicator Spatial Autocorrelaetion
dev.off()


#moran's I plot
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/MoranI_Plot.png", width = 800, height=600)
moran.plot(income.pm25$Income, vancouver.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Values at Location i", 
           ylab="Values at Neighbourhood i", quiet=NULL)
dev.off()


#---------------------------------------------------------------- Spatial Interpolation: Ord. Kriging ----------------------

# #Ordinary Kriging
# 
# f.0 <- as.formula(PM25AGG ~ 1) 
# 
# var.smpl <- variogram(f.0, spSample, cloud = FALSE) #, cutoff=10000, width=89900)
# 
# dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
#                           vgm(model="Exp"))
#                           #vgm(psill=1.9e-05, model="Gau", range=17700, nugget=0))
# 
# #here we actually define how the semivariogram (line) behaves, compared to our data
# plot(var.smpl, dat.fit)
# 
# # Perform the krige interpolation (note the use of the variogram model)
# dat.krg <- krige(f.0, spSample, grd, dat.fit)
# 
# # Convert kriged surface to a raster object for clipping
# r <- raster(dat.krg)
# r.m <- mask(r, income.pm25)
# 
# # Plot the krige surface
# #tmap_mode("view")
# 
# tm_shape(r.m) + 
#   tm_raster(n=10, palette="Oranges", alpha = 1, 
#             title="Predicted PM2.5 \n(in ppm)") +
#   tm_shape(spSample) + 
#   tm_dots(col = "PM25AGG", size=0.02) +
#   tm_legend(legend.outside=TRUE)
# 
# # Plot the Variance surface
# r   <- raster(dat.krg, layer="var1.var")
# r.m <- mask(r, income.pm25)
# 
# tm_shape(r.m) + 
#   tm_raster(n=7, palette ="Reds",
#             title="Variance map \n(ppm\u00B2)") +tm_shape(spSample) + tm_dots(size=0.02) +
#   tm_legend(legend.outside=TRUE)
# 
# # Plot the 95% CI map
# r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
# r.m <- mask(r, income.pm25)
# 
# tm_shape(r.m) + 
#   tm_raster(n=7, palette ="Reds",
#             title="95% CI map \n(in ppm)") +tm_shape(spSample) + tm_dots(size=0.02) +
#   tm_legend(legend.outside=TRUE)


##--------------------------------------------------------------- Spatial Interpolation: with Polynomial Trends


#------------------------------------------- Define the 2nd order polynomial equation
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
spSample$X <- coordinates(spSample)[,1]
spSample$Y <- coordinates(spSample)[,2]

# Run the regression model
lm.2 <- lm(f.2, data=spSample)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

tmap_mode("plot")

# Clip the interpolated raster to Vancouver
r   <- raster(dat.2nd)
r.m <- mask(r, income.pm25)

# Plot the 2nd-order polynomial map
tm_shape(r.m) + 
  tm_raster(n=10, palette="Oranges",
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(spSample) + tm_dots(size=0.02) +
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))



#------------------------------------------------------------------------- Spatial Interpolation: Kriging (run LAST) 

#-------------------------------------------------- Universal Kriging -------------------------------

f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

var.smpl <- variogram(f.2, spSample, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model="Sph"))

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/SecondOrderSemivariogram.png", width = 800, height=600)
plot(var.smpl, dat.fit)
dev.off()

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg2 <- krige(f.2, spSample, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg2)
r.m_Krig <- mask(r, income.pm25)

# Plot the 2nd-order Kriged map

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/SecondOrderKrigingSurface_8values(2).png", width = 800, height=600)
tm_shape(r.m_Krig) + 
  tm_raster(n=8, 
            palette="-RdBu", 
            midpoint = 0,
            alpha = 1, 
            title="Predicted PM2.5 (ppm)") +
  tm_shape(spSample) + 
  tm_dots(title = "PM2.5 data points", col = "PM25AGG", palette = "Reds", size=0.2)+
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))
dev.off()

# Plot the 2nd-Order Variance Surface
r   <- raster(dat.krg2, layer="var1.var")
r.m <- mask(r, income.pm25)

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/SecondOrderVarianceSurface.png", width = 800, height=600)
tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds", alpha = 1,
            title="Variance map \n(ppm\u00B2)") +tm_shape(spSample) + tm_dots(size=0.02) +
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("right","top"))
dev.off()

# Plot the 2nd-Order 95% CI surface
r   <- sqrt(raster(dat.krg2, layer="var1.var")) * 1.96
r.m <- mask(r, income.pm25)
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/SecondOrderCISurface.png", width = 800, height=600)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map \n(in ppm)") +tm_shape(spSample) + tm_dots(size=0.02) +
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("right","top"))
dev.off()


#--------------------------- Combine outputs from Spatial Interpolation with Income Data.

#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the max.
#essentially choosing the max PM2.5 value from each Dissemination Area to represent the DA.   


#If you have too many cells, you can reduce the number by aggregating values
step.1 <- aggregate(r.m_Krig, fact=1, fun=mean)
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/Agg.png", width = 800, height=600)
plot(step.1)
dev.off()

#if you didn't want to optimize the above (i.e.increase your cell size by (fact)^2), use the raster you created from your interpolated surface. 
#Convert the raster dataset to points
step.2 <-  rasterToPoints(r.m_Krig,fun=NULL, spatial=FALSE, crs=crs(pm25.spatial))
step.2 <- as.data.frame(step.2) #convert the point dataset to a spatial dataframe
# colnames(step.2)[3] <- "PM25"
Coords <- step.2[,c("x", "y")]  #assign coordinates to a new object
crs <- crs(census.tracts) #utilize an existing projection
step.3 <- SpatialPointsDataFrame(coords = Coords, data = step.2, proj4string = crs) #create a spatial points dataframe
step.4 <- aggregate(x=step.3,by=income.tracts, FUN=mean) #aggregate points into census tracts
step.5 <- intersect(step.4,income.tracts)  #get the intersection of step.4 with the income.tracts dataset (this will take a while) 

pm.income.poly <- step.5
#----------------------------------------------------------- Linear Regression


png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/PM25RegressionMap.png", width = 800, height=600)
tm_shape(pm.income.poly)+
  tm_polygons(col="var1.pred", palette="-RdBu", lwd = 0.5,
              title="Predicted PM2.5 \n(in ppm)") +
    tm_legend(legend.outside=TRUE)

#Let's say your dataset with both PM2.5 and Income are stored in a dataset called pm.income.poly.

#Plot income and PM2.5 from the pm.income.poly dataset you created
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/PM25RegressionMap_WITHZEROS.png", width = 800, height=600)
plot(pm.income.poly$Income~pm.income.poly$var1.pred) #There are a lot of 0's in this dataset. If you decide to remove them, use the following line:
dev.off()

new_DF <- pm.income.poly #saving the original dataframe 
NoNa.pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$var1.pred),] #removing NA values
NoNa.pm.income.poly <- subset(NoNa.pm.income.poly, var1.pred >= minPM25 & var1.pred <= maxPM25) #taking values that are smaller than the max and greater than the minimum PM2.5 value
#pm.income.poly <-  pm.income.poly[pm.income.poly$PM25 != 0, ]
names(NoNa.pm.income.poly)[3] <- "PM25" #changing column names because who knows what var1.pred means

#Now plot the data again
png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/PM25RegressionMap_NoNAorZEROS.png", width = 800, height=600)
plot(NoNa.pm.income.poly$Income~NoNa.pm.income.poly$PM25)
dev.off()

#Perform a linear regression on the two variables. You should decide which one is dependent.
NoNa.lm.model <- lm(NoNa.pm.income.poly$PM25~NoNa.pm.income.poly$Income)
#lm.model <- lm(pm.income.poly$Income~pm.income.poly$PM25)


#Add the regression model to the plot you created
abline(NoNa.lm.model)
#Get the summary of the results
summary(NoNa.lm.model)

#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(NoNa.lm.model))

#Then add the residuals to your spatialpolygon dataframe
NoNa.pm.income.poly$residuals <- residuals.lm(NoNa.lm.model)

#Observe the result to make sure it looks correct
head(NoNa.pm.income.poly)


# #Now, create choropleth map of residuals
# resids <- NoNa.pm.income.poly$residuals
# shades <- auto.shading(resids, n=6, cols = brewer.pal(6, 'Greens'))
# choropleth(income.tracts, resids, shades) #map the data with associated colours
# choro.legend(-123.5075, 49.15346, shades) #add a legend (you might need to change the location


#--------------------------------------------------------------------------------------------------- Alternative Mapping Technique

#WHEN MAPPING RESIDUALS: DO THE NEGATIVES CANCEL OUT THE POSITIVES? 

map_PM25Res <- tm_shape(NoNa.pm.income.poly) + 
  tm_polygons(col = "residuals", 
              title = "PM2.5 Regression Residuals\nVancouver, British Columbia", 
              style = "jenks", #suggested by R for larger datasets
              palette = "Greens", n = 6, lwd = 0.5,
              midpoint = NA)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("left","bottom"))


png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/PM25ResidualsMapNEW.png", width = 800, height=600)
map_PM25Res #Local Indicator Spatial Autocorrelaetion
dev.off()


#---------------------------------------------------------------------------------------------------- Moran's I on residuals


#Weighted Matrix, spatially based on Census Tract information

van.res.nb <- poly2nb(NoNa.pm.income.poly) #using study area (defined by income) tracts, queens weighting -> defining a polygon's neighbourhood, up to eight surrounding points
van.res.net <- nb2lines(van.res.nb,coords=coordinates(NoNa.pm.income.poly))

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/ResidualsQueensNeighbourhood.png", width = 800, height=600)
tm_shape(NoNa.pm.income.poly) +
  tm_borders(col='lightgrey') +
  tm_shape(van.res.net) +
  tm_lines(col='blue') #creating a visualization of the Queens case weighted mean
dev.off()


#creating the weight matrix (i.e a bunch of 1's and 0's, defining if they're neighbours or not)
van.res.lw <- nb2listw(van.res.nb, zero.policy = TRUE, style = "W")
print.listw(van.res.lw, zero.policy = TRUE)


#-------------------------------------------------------------------------------------------------------------------Global Moran's I

mi <- moran.test(NoNa.pm.income.poly$residuals, van.res.lw, zero.policy = TRUE) #based off vancouver.lw weighted matrix
mi

moran.range <- function(lw) { 
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(van.res.lw) #we are going to need these values to calculate the Z-score


mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- ((mI-eI)/(sqrt(var)))
z


#-------------------------------------------------------------------------------------------------------------------- Local Moran's I
residual_lisa.test <- localmoran(NoNa.pm.income.poly$residuals, van.res.lw) #based off va.res.lw weighted matrix

#creating columns for expected values
NoNa.pm.income.poly$Ii <- residual_lisa.test[,1] #this describes the expected values we are going to get: I, E,Var, Z,P
NoNa.pm.income.poly$E.Ii<- residual_lisa.test[,2]
NoNa.pm.income.poly$Var.Ii<- residual_lisa.test[,3]
NoNa.pm.income.poly$Z.Ii<- residual_lisa.test[,4]
NoNa.pm.income.poly$P<- residual_lisa.test[,5]

#-------------------------------------------------------------------------------------------------------------------- mapping the local I

map_res.LISA <- tm_shape(NoNa.pm.income.poly) +
  tm_polygons(col = "Ii",
              title = "Local Moran's I value of Median Income,\nFisher Classification",
              style = "fisher", #suggested by R for larger datasets
              palette = "-RdBu", n = 5, lwd = 0.5,
              midpoint = NA)

#tmap_mode("view")
tmap_mode("plot")

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/res.LISA_MoransIValue_IncomeMap.png", width = 800, height=600)
map_res.LISA #Local Indicator Spatial Autocorrelaetion
dev.off()

#plot based on p - values : which indicate clustering (low p value) - define breaks
map_res.LISA_p <- tm_shape(NoNa.pm.income.poly) +
  tm_polygons(col = "P",
              title = "Local Moran's I p-value of Median Income,\nFisher Classification",
              style = "fisher", #suggested by R for larger datasets
              palette = "-RdBu", n = 5, lwd = 0.5)

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/resLISA_PValue_IncomeMap.png", width = 800, height=600)
map_res.LISA_p #Local Indicator Spatial Autocorrelaetion
dev.off()


#----------------------------------------------------------------------- GWR

####Geographically Weighted Regression

#Let's say you are continuing with your data from the regression analysis. 

#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(NoNa.pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
NoNa.pm.income.poly$X <- pm.income.poly.coords[,1]
NoNa.pm.income.poly$Y <- pm.income.poly.coords[,2]
head(NoNa.pm.income.poly)

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(NoNa.pm.income.poly$PM25~NoNa.pm.income.poly$Income, 
                        data=NoNa.pm.income.poly, coords=cbind(NoNa.pm.income.poly$X,NoNa.pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(NoNa.pm.income.poly$PM25~NoNa.pm.income.poly$Income, 
                data=NoNa.pm.income.poly, coords=cbind(NoNa.pm.income.poly$X,NoNa.pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
NoNa.pm.income.poly$localr <- results$localR2

NoNa.pm.income.poly2 <- subset(NoNa.pm.income.poly,localr>=0)

# shades <- auto.shading(local.r.square, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, local.r.square, shades) #map the data with associated colours
# choro.legend(-123.5075, 49.21703, shades) #add a legend (you might need to change the location)

LocalR2 <- tm_shape(NoNa.pm.income.poly2) + 
  tm_polygons(col = "localr", 
              title = "Local R\u00B2 of PM2.5\nVancouver, British Columbia", 
              style = "jenks",
              palette = "-Greys", n = 3, alpha = 0.7, lwd = 0.5, border.alpha = 0.25,
              midpoint = NA)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("right","TOP"))+
  tm_legend(legend.outside=TRUE)

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/LocalR2MapNONA.png", width = 800, height=600)
LocalR2 #Local Indicator Spatial Autocorrelaetion
dev.off()


#Time for more magic. Let's map the coefficients
NoNa.pm.income.poly$coeff <- results$NoNa.pm.income.poly.Income

# #Create choropleth map of the coefficients
# local.coefficient <- NoNa.pm.income.poly$coeff
# shades <- auto.shading(local.coefficient, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, local.coefficient, shades, cex = 0.8) #map the data with associated colours
# choro.legend(-123.5075, 49.21703, shades) #add a legend (you might need to change the location)

GWR_coeff <- tm_shape(NoNa.pm.income.poly) +
  tm_polygons(col = "coeff",
              title = "GWR Coefficients Map\nVancouver, British Columbia",
              style = "jenks", 
              palette = "-RdBu", n = 5, lwd = 0.5, border.alpha = 0.5,
              midpoint = NA)+
  tm_scale_bar(width = 0.15,text.size = 1.5, position = c("right","TOP"))+
  tm_legend(legend.outside=TRUE)

png(filename = "D:/GEOG418/FinalProject/geog418-518-2019-finalproject-master/output/GWRCoeffMap.png", width = 800, height=600)
GWR_coeff
dev.off()
#-------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------- Point Pattern Analysis 

#------- clean dat data

View(spSample)

spSample <- spTransform(spSample, CRS("+init=epsg:3005")) #projecting data in meters

spSample$X0 <- coordinates(spSample)[,1]
spSample$Y0 <- coordinates(spSample)[,2]

#create an "extent" object which can be used to create the observation window for spatstat
spSample.ext <- as.matrix(extent(spSample)) 

#observation window
window <- as.owin(list(xrange = spSample.ext[1,], yrange = spSample.ext[2,]))

#create ppp oject from spatstat
spSample.ppp <- ppp(x = spSample$X0, y = spSample$Y0, window = window)

#--------------------------- Kernel Density Estimator 


kde.100 <- density(spSample.ppp, sigma = 2040, at = "pixels", eps = c(100, 100))
kde.SG <- as(kde.100, "SpatialGridDataFrame")

  
names(kde.SG) <- c("sigma100"," sigma250")
#plot
spplot(kde.SG)

#can see how the bandwidth selection influences the density estimates
summary(kde.SG)

#use cross-validation to get the bandwidth that minimizes MSE
bw.d <- bw.diggle(spSample.ppp)
#plot the "optimal" bandwidth
plot(bw.d, ylim=c(-10, 10), main= "Cross-Validation KDE Bandwidth Optimization Tool")

#density using the cross-validation bandwidth
kde.bwo <- density(spSample.ppp, sigma = bw.d, at = "pixels", eps = c(50, 50))
plot(kde.bwo)



#ENJOY! Thanks for the great term,
# Cole Glover / V00838832 
