library("ggplot2", lib.loc="~/R/win-library/3.2")
library("ggmap", lib.loc="~/R/win-library/3.2")
library("akima", lib.loc="~/R/win-library/3.2")
library("rgdal", lib.loc="~/R/win-library/3.2")
library("rgeos", lib.loc="~/R/win-library/3.2")
library("maptools", lib.loc="~/R/win-library/3.2")
source('~/R Projects/UmatillaReds/R/spatialUtils.R')


dsn <- "C:/Users/goff/Documents/R Projects/HouseMap/shapefiles"
points <- readOGR(dsn=dsn, layer="Point_ge")
coordnames(points) = c("long", "lat")
names(points)[names(points) == "GNSS_Heigh"] = "z"

# replace some bad data
points[30, "Comment"] = "elevation"
points[5, "z"] = 1523.836

#read in the edges file (road and trail)
edges = readOGR(dsn = dsn, layer = "Line_gen")
coordnames(edges) = c("long", "lat")

dap = distAlongSpatialLines(points, edges, withAttrs = F)
skidMarks = lapply(
  row.names(dap), 
  function(.i) {
    Lines(
      list(
        Line(
          coordinates(
            rbind(
              geometry(points[.i,]),
              geometry(dap[.i,])
            )
          )
        )
      ),
      .i
    )
  }
)


skidMarks = SpatialLines(skidMarks, proj4string = points@proj4string)

plot(points)
lines(edges)
points(dap, col = "red")
lines(skidMarks, col = "red")

head(dap, 10)
