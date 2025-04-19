
                                            
#### ALS-BASED FOREST INVENTORY WITH LIDR ####



#==================================
# Initialize R
#==================================

# Clear memory of old stuff
rm(list = ls())

# install.packages("sf")
# install.packages("lidR")
# install.packages("terra")
#install.packages("microbenchmark")

# Load installed libraries into active use
require(sf)
require(terra)
require(lidR)
require(microbenchmark)

# Check if the help page of the following command shows up. 
# If yes, the installation is ok.
#polygon_metrics

# Set working directory path - all your files should be here
setwd("C:\\Users\\Hp\\OneDrive - University of Eastern Finland\\Desktop\\GIT-HUB\\Als-based-forest-inventory-with-lidar")

# Load a set of custom functions
source("function.R")


#==================================
# Read data sets into R
#==================================

# Read background and plot aerial image
ortho <- rast("kiihtelys_ortho.tif")
plotRGB(ortho)

# Read field plot shapefile from the working directory
field <- st_read("ars_plots.shp")

# Show field plot locations over the aerial image
plot(field, add=TRUE, col="red")

# View the attributes of the field plots
View(field)

# Read laser data file into R's LAS object
las <- lidR::readLAS("kiihtelys.laz")

# View a summary of the LAS object
summary(las)

# Show the records of the LAS object
tail(las@data)

# View the a table of classifications of the las object
table(las$Classification)


# Add echo type to User Data, overwrite the previous las object
las <- add_echo_type(las)

# View the updated User Data field in tabular form
table(las$UserData)

# Check laser data integrity
las_check(las)

#========================================
# Visualization 
#========================================

# Clip a smaller portion of the las file for a better view
subarea <- clip_rectangle(las, 665000, 6934700, 665300, 6935000)

# Plot the smaller area in different ways
plot(subarea, legend=T) # Elevation
#plot(subarea, color="UserData", legend=T) # User Data
#plot(subarea, color="Classification", legend=T) # Classification
#plot(subarea, color="Intensity", legend=T) # Intensity

# Re-plot Intensity with rescaled colors
plot(subarea, color="Intensity", breaks = "quantile", nbreaks = 50, legend=T)

# Clip point clouds of plots from the LAS object
 cliped <- clip_roi(las, field)

# View point clouds of plots
#plot(cliped[[10]]) #plot 1 point clouds


#=================================================
# Height normalization and removal of outliers
#=================================================

# Construct a Digital Terrain Model
# Set output file, input file, spatial resolution, and interpolation method
dtm <- rasterize_terrain(las, res = 1, tin())


# View the DTM in 2D and 3D hillshade
plot(dtm)
plot_dtm3d(dtm)

# Height normalize LAS object by subtracting the DTM
dzlas <- las - dtm

# Remove the non-normalized LAS object and the DTM to free memory
rm(las)
rm(dtm)

# Reset graphical parameters
dev.off()

# View normalized Z and intensity histograms
hist(dzlas$Z)
hist(dzlas$Intensity)

# Remove noise echoes
dzlas <- classify_noise(dzlas, ivf(res = 5, n = 2))
dzlas <- filter_poi(dzlas, Classification != LASNOISE)

# Remove echoes with unrealistic intensities
dzlas <- filter_poi(dzlas, Intensity < 300)

# Review echo histograms
hist(dzlas$Z)
hist(dzlas$Intensity)

# Set negative elevations as zero
dzlas$Z[dzlas$Z < 0] <- 0

max(dzlas$Z) # View Z maximum
min(dzlas$Z) # View Z minimum

max(dzlas$Intensity) # View Intensity maximum
max(dzlas$Intensity) # View Intensity minimum

#============================================
# Computation of predictor variables
#============================================

# Compute preditor variables for each plot into data frame d
# Inputs: LAS file, function that computes the variables, plot polygons
d <- polygon_metrics(dzlas, ~uef_metrics(Z, Intensity, UserData), field)



tail(d)      # View structure of the predictor file
class(d)     # Note that it can be treated as a georeferenced data frame

# Merge with field plots using column bind
d <- cbind(field, d)

# Convert into a normal data frame
d <- st_set_geometry(d, NULL)


tail(d)      # View structure of the predictor file
class(d)     # Note that it can be treated as a georeferenced data frame
names(d)
attach(d)

# Make a scatter plot of max echo height (f_hmax) and dominant height (hdom)
# If the predictors are OK, there should be a high correlation
plot(f_hmax, hdom, xlim=c(0,30), ylim=c(0,30)); abline(0,1)
cor(f_hmax, hdom)

#=================================================
# Construct a regression model for plot volume
#=================================================

# Make your own plot level volume model 
summary(m <- lm(v~f_hmean + f_hstd, data=d))
plot(fitted(m), resid(m) ,xlab= "Fitted values", ylab= "Residuals", abline(h=0))
rmse(d$v, fitted(m))
# exh_var_search(d, 4, "v", 12, 68) # variable search algorithms:
#V1     V2    V3    V4  RMSE
#160   f_hmean f_hstd f_h80 f_h70 0.072
#26148 f_hmean l_hstd l_p80 l_p40 0.072
#1453  f_hmean f_hstd l_p80 l_p40 0.073
#27662 f_hmean  l_p80 l_p20 l_p10 0.073
#26149 f_hmean l_hstd l_p80 l_p30 0.073


summary(m.1 <- lm(v ~ f_hmean +  f_h80 +  f_h70))
summary(m.2 <- lm(v ~ f_hmean +   l_p80 + l_p40))
summary(m.3 <- lm(v ~ f_hmean + f_hstd + l_p80 + l_p40))
summary(m.4 <- lm(v ~ f_hmean + l_hstd + l_p80 +  l_p30))
summary(m.1)
rmse(d$v, fitted(m.1))
plot(fitted(m.2), resid(m.2) ,xlab= "Fitted values", ylab= "Residuals", abline(h=0))
plot(predict(m.2), d$v, xlim=c(0,400), ylim=c(0,400),xlab="Predicted", ylab="Observed", abline(0,1))
bias(d$v, fitted(m.2))
rmse(d$v, fitted(m.2))
relrmse (d$v, fitted(m.2))
relbias(d$v, fitted(m.2))


#Transformation:
#d$f_mean.1 <- sqrt(d$f_hmean)
#d$f_h80.1 <- sqrt(d$f_h80)
#d$f_h70.1 <- sqrt(d$f_h70)

#transformed model:
m.t <- lm(v ~ f_mean.1 + f_h80.1 + f_h70.1, data = d)
summary(m.t)
rmse(d$v, fitted(m.t))
plot(fitted(m.t), resid(m.t) ,xlab= "Fitted values", ylab= "Residuals", abline(h=0))
plot(fitted(m.t), resid(m.t) ,xlab= "Fitted values", ylab= "Residuals", abline(h=0))
plot(predict(m.t), d$v, xlim=c(0,400), ylim=c(0,400),xlab="Predicted", ylab="Observed", abline(0,1))


#==================================================
# Wall-to-wall prediction
#==================================================

# Compute a raster stack of ALS metrics for the whole area
# Parameters: LAS object, function to compute, resolution
abagrid <- pixel_metrics(dzlas, ~uef_metrics(Z,  Intensity, UserData), res = 15)

# Plot the f_hmean from the raster stack
plot(abagrid$f_hmean)

# Predict volume using the raster stack and your model object
 
# Save as v_pred
v_pred <- 10557.025 + 27.758 * abagrid$f_hmean   -15344.057 * abagrid$l_p80 + 4773.840* abagrid$l_p40
predict_values<- predict(abagrid, m.2)
# Force areas where top echo is height is < 2 m as 0
v_pred[abagrid$f_hmax < 2] <- 0

# Reset graphics and plot volume grid
#dev.off()
plot(v_pred, main="Volume m3/ha")
summary(v_pred)
v_pred[v_pred < 0] <- 0 #re-plot it after the correction to negative values
writeRaster(v_pred, 'Volume_pred.tif', overwrite= T)
# Save volume grid as a .tif file and view it in a GIS
terra::writeRaster(v_pred, "pred.v.1.tif", overwrite=T)

# Free memory
rm(abagrid)

#============================================
# Construction of a canopy height model
#============================================

# Interpolate a CHM
# Parameters: LAS object, spatial resolution, interpolation algorithm
chm <- rasterize_canopy(dzlas, res = 0.5, algorithm = pitfree())

# Plot CHM
plot(chm)

# Add field plot polygons 
plot(field, add=T, col="black")

# Crop a smaller area from the CHM for closer inspection and plot it
subchm <-crop(chm, ext(664030, 664090, 6935380, 6935440))
plot(subchm)


#============================================
# Individual tree detection
#============================================

# Define a function for computing a height-dependent moving window size
window_size <- function(height){ 1.2 + 0.003*height^2 } 

# Compute a georeferenced data frame with tree top coordinates and heights
tops <- locate_trees(chm, lmf(window_size))

tail(tops) # View output


# Visualize tree locations
plot(subchm) # Re-plot the small-area CHM

# Re-add plot polygons as transparent
plot(field, add=T, lwd=2, col=rgb(0,0,0,0)) 

# Add detected tree locations on top
plot(tops, add=T)


# Apply height correction model with the tree top data frame
tops$h <- 1.2 + tops$Z*1.01

tail(tops) # View the result

#============================================
# Construction of a tree level volume model
#============================================

dev.off() # Reset graphics

# Read in tree file
trees <- read.csv("ars_trees.csv", header=T, sep=",", dec=".")
dim(trees)


# Make your own tree level volume model 
model.t <- lm(v ~ h , data = trees)
summary(model.t)
rmse(trees$v, fitted(model.t))
plot(fitted(model.t), resid(model.t) ,xlab= "Fitted values", ylab= "Residuals", abline(0,0))
plot(predict(model.t), trees$v,xlab="Predicted", ylab="Observed")
v_pred <- predict(model.t)
plot(trees$h.t, trees$v)
x <- seq(0,30,1) # Creates numeric sequence 0, 1, 2, …, 30
y <- exp(-9.61 + 1.92 * x + (0.0242)^2/2)


#transform model: log transformation:
trees$h.t <- log(trees$h)
trees$log.v<- log(trees$v)
model.t.t <- lm(log.v ~  h.t, data = trees)
plot(fitted(model.t.t), resid(model.t.t),abline(0,0))
summary(model.t.t)
rmse(trees$log.v, fitted(model.t.t))
log.pred<-predict(model.t.t)
#back-transformation:
pred_v <- exp(log.pred+ (0.4443)^2/2)

#Scatter plot of height and volume:

### One question:
plot(  trees$h, trees$v, xlim = c (0,30), ylim = c(0,3), xlab = "Trees Volume", ylab = "Trees height")
x <- seq(0,30,1) # Creates numeric sequence 0, 1, 2, …, 30
y <- exp(-11.56 + 3.52 * log(x) + (0.4443)^2/2)
  #(trees$d13) + 1.170e-04 * (trees$h)
 # Apply your model coefficients to the sequence x
lines(x, y, col="red")

# Apply your model to predict a new column in the tree top file
rmse(trees$v, fitted(model.t.t))
plot(fitted(model.t.t), resid(model.t.t) ,xlab= "Fitted values", ylab= "Residuals", abline(h=0))
plot(predict(model.t.t), resid(model.t.t), xlim=c(0,400), ylim=c(0,400),xlab="Predicted", ylab="Observed", abline(0,1))

#Adding volumes:
tops$v_itd <- exp(-11.56 + 3.52 * log(tops$h) + ((0.4443)^2/2))


#============================================
# Plot level accuracy assessment for ITD
#============================================

# Convert the detected trees to a different spatial format
topsf <- st_as_sf(tops)         
 #class(tops) 
# Set the coordinate system of this object to be the same as the field plots
st_crs(tops) <- 3067 

# Make a spatial overlay of field plots and detected trees 
plot_tops <- st_join(field,tops) 

# View the last rows of the new object with plot-labeled itd trees
head(plot_tops)

# Compute a new column called v_itd in your FIELD PLOT TABLE using function 
# tapply with your plot-labeled itd trees. Tapply applies a function (here sum())
# to the values of a given column (v_itd) using a given grouping variable (plot).
field$v_itd <- tapply(plot_tops$v_itd, plot_tops$plot_id, sum)

# Repeat, but now compute the plot level stem number by using the length() function
# instead of sum() and save it as n_itd
field$n_itd <- tapply(plot_tops$n, plot_tops$plot_id, length)

# Scale v_itd and n_itd to per hectare using plot radius
field$v_itd <- field$v_itd * 10000
field$n_itd <- field$n_itd * 10000

# View the last rows of the field plot table
tail(field)

# Make scatter plots of plot level predicted vs observed volume and stem number
plot(field$v_itd, field$v, xlim = c(0, 475), ylim = c(0, 475))

# Compute plot level RMSE-% and bias-% of both volume and stem number
relrmse(field$v_itd, field$v)
relbias(field$n_itd, field$n)

