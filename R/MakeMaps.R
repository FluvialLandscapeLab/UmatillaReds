library("ggplot2", lib.loc="~/R/win-library/3.2")
library("ggmap", lib.loc="~/R/win-library/3.2")
library("akima", lib.loc="~/R/win-library/3.2")
library("rgdal", lib.loc="~/R/win-library/3.2")
library("rgeos", lib.loc="~/R/win-library/3.2")
library("maptools", lib.loc="~/R/win-library/3.2")
source('~/R Projects/UmatillaReds/R/spatialUtils.R')

sampPtSpacing = 100

#read in the points file and replace a couple of names...
dsn <- "C:\\Users\\goff\\Documents\\R Projects\\UmatillaReds\\SpatialData\\FloodplainBoundary"
FPN <- readOGR(dsn=dsn, layer="UmatillaFloodplainBoudary_North")
FPS <- readOGR(dsn=dsn, layer="UmatillaFloodplainBoudary_South")
dsn <- "C:\\Users\\goff\\Documents\\R Projects\\UmatillaReds\\SpatialData\\ValleyCenterLine"
FPC <- readOGR(dsn=dsn, layer="ValleyCenterLine_Dissolved")

#FPC is not a single line.  It has a couple of spurs.  Fix this.
FPC = simplifySpatialLines(FPC, minSpurLength = 50)

plot(FPN)
lines(FPS)    
lines(FPC, col = "red")

samplingPoints = repeatAlongLines(FPC@lines[[1]], sampPtSpacing, 0, proj4string = FPC@proj4string)
points(samplingPoints)

samplingPoints@data[,"FPWidth"] = as.vector(gDistance(FPN, samplingPoints, byid = T) + gDistance(FPS, samplingPoints, byid = T))

head(samplingPoints)

plot(samplingPoints@data$distance, samplingPoints@data$FPWidth, type="l")
samplingPoints@data
