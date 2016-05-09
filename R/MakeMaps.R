library("ggplot2", lib.loc="~/R/win-library/3.2")
library("ggmap", lib.loc="~/R/win-library/3.2")
library("akima", lib.loc="~/R/win-library/3.2")
library("rgdal", lib.loc="~/R/win-library/3.2")
library("rgeos", lib.loc="~/R/win-library/3.2")
library("maptools", lib.loc="~/R/win-library/3.2")

source(paste0(getwd(),'/R/spatialUtils.R'))
source(paste0(getwd(),'/R/reddAndValleyWidthAnalysis_Utils.R'))

sampPtSpacing = 0.05

dsn <- "C:\\Users\\goff\\Documents\\R Projects\\UmatillaReds\\SpatialData\\FloodplainBoundary"
FPN <- readOGR(dsn=dsn, layer="UmatillaFloodplainBoudary_North")
FPS <- readOGR(dsn=dsn, layer="UmatillaFloodplainBoudary_South")

dsn <- "C:\\Users\\goff\\Documents\\R Projects\\UmatillaReds\\SpatialData\\ValleyCenterLine"
FPC <- readOGR(dsn=dsn, layer="ValleyCenterLine_Dissolved")
#FPC is not a single line.  It has a couple of spurs.  Fix this.
FPC = simplifySpatialLines(FPC, minSpurLength = 50)

dsn <- "C:\\Users\\goff\\Documents\\R Projects\\UmatillaReds\\SpatialData\\Redds"
redds <- readOGR(dsn = dsn, layer = "redds_mainstemUmatilla")

#convert projection to km
newproj4string = CRS("+proj=utm +zone=11 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")
targSpNames = c("FPN", "FPS", "FPC", "redds")
for(tName in targSpNames) {
  assign(tName, spTransform(get(tName), newproj4string))
}

redds@data = redds@data[c("ANADROMOUS", "GPS_DATE")]
redds = distAlongSpatialLines(redds, FPC)

samplingPoints = repeatAlongLines(FPC@lines[[1]], sampPtSpacing, 0, proj4string = FPC@proj4string)
samplingPoints@data[,"FPWidth"] = as.vector(gDistance(FPN, samplingPoints, byid = T) + gDistance(FPS, samplingPoints, byid = T))

tail(samplingPoints)

# plot(FPN)
# lines(FPS)    
# lines(FPC, col = "blue", lwd = 1)
# points(samplingPoints, col = rainbow(100))
# points(redds)

##### Parameters ################

# downstream point on floodplain to start including data
startKm = 110

samplingPoints = samplingPoints[samplingPoints@data$distance >= startKm,]

boundsListList = sapply(
  c(down = -1, up = 1), 
  segBounds, 
  location = samplingPoints$distance, 
  distance = samplingPoints$FPWidth,
  simplify = F
)

plot(samplingPoints@data$distance, samplingPoints@data$FPWidth, ty = "l")
troughIndexes = which(calcPeaks(samplingPoints@data$FPWidth, fun = min, nPreAndPost = 4))
points(samplingPoints@data$distance[troughIndexes], samplingPoints@data$FPWidth[troughIndexes])
peakIndexes = which(calcPeaks(samplingPoints@data$FPWidth, nPreAndPost = 4))
points(samplingPoints@data$distance[peakIndexes], samplingPoints@data$FPWidth[peakIndexes])

downstreamPeakIdx = sapply(troughIndexes, function(.i) peakIndexes[max(which(peakIndexes<.i))])
upstreamPeakIdx = sapply(troughIndexes, function(.i) peakIndexes[min(which(peakIndexes>.i))])

zoneBounds = data.frame(start = downstreamPeakIdx, trough = troughIndexes, end = upstreamPeakIdx, row.names = NULL)
# duplicates in start and end indices indicate a location where a minor trough
# separates narrowing from widening zone.  Delete these.
zoneBounds = zoneBounds[!(zoneBounds$start %in% zoneBounds$start[duplicated(zoneBounds$start)]),]
# get rid of NAs
zoneBounds = zoneBounds[!apply(zoneBounds, 1, function(x) any(is.na(x))), ]

zoneBounds["11", "end"] = 141

zoneWidths = apply(zoneBounds, 1, function(x) list(wide = samplingPoints$FPWidth[x[1]:x[2]], narrow = samplingPoints$FPWidth[x[2]:x[3]]))

zoneWidths = lapply(
  zoneWidths, 
  function(.l) {
    targLen = min(length(.l$narrow),length(.l$wide))
    return(lapply(.l, "[", 1:targLen))
  }
)


fracNarrowing = sapply(zoneWidths, function(.l) (max(.l$narrow) - min(.l$narrow)) / max(.l$narrow))
keepEm = fracNarrowing > 0.33
zoneWidths = zoneWidths[keepEm]
zoneBounds = zoneBounds[keepEm,]


newBounds = sapply(
  1:nrow(zoneBounds),
  function(.i) {
    c(
      start = zoneBounds[.i,"trough"] - (length(zoneWidths[[.i]]$narrow)-1),
      end = zoneBounds[.i,"trough"] + (length(zoneWidths[[.i]]$wide)-1)
    )
  }
)

### Useful to make sure newBounds is correct...
# cbind(zoneBounds$start, zoneBounds$end, newBounds["start",], newBounds["end",])

zoneBounds$start = newBounds["start",]
zoneBounds$end = newBounds["end",]


zoneIndexes = mapply(":", zoneBounds[,"start"], zoneBounds[,"trough"], SIMPLIFY = F)
zoneIndexes = c(zoneIndexes, mapply(":", zoneBounds[,"trough"], zoneBounds[,"end"], SIMPLIFY = F))

reddCounts = t(apply(zoneBounds, 1, localReddCount, redds=redds, samplingPoints = samplingPoints))

t.test(reddCounts[,"wide"], reddCounts[,"narrow"], paired = T)

color = rep(c("blue", "red"), each = nrow(zoneBounds))

plot(samplingPoints$distance, samplingPoints$FPWidth, col = "grey", lwd = 2, ty = "l")
dummy = lapply(
  1:length(zoneIndexes),
  function(.i) {
    lines(
      x = samplingPoints$distance[zoneIndexes[[.i]]],
      y = samplingPoints$FPWidth[zoneIndexes[[.i]]],
      lwd = "2",
      col = color[.i]
    )
  }
)

plot(samplingPoints[250:350,], col = "white")
lines(FPN)
lines(FPS)
dummy = lapply(
  1:length(zoneIndexes),
  function(.i) {
    points(
      samplingPoints[zoneIndexes[[.i]],],
      col = color[.i]
    )
  }
)
points(redds, pch = ".")

