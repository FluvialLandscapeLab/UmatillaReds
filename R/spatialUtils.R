.fractionBetween = function(xy, bracket) {
  frac = (xy-bracket[1,])/(bracket[2,]-bracket[1,])
  return(frac[!is.nan(frac)][1])
}

spatialLinesFromEndPoints = function(coords1, coords2, proj4string) {
  if(!is.matrix(coords1) || !is.matrix(coords2)) stop("Coordinates must be of type 'matrix'")
  if(nrow(coords1) != nrow(coords2)) stop("Coordinate matricies must have the same number of rows.")
  if(ncol(coords1) != ncol(coords2) || ncol(coords1) != 2) stop("Coordinate matricies must have 2 columns each.")
  segList = lapply(
    1:nrow(coords1), 
    function(.i) {
      Lines(
        list(
          Line(
            rbind(
              coords1[.i,],
              coords2[.i,]
            )
          )
        ),
        .i
      )
    }
  )
  return(SpatialLines(segList, proj4string = proj4string))
}

contourLinesSGDF = function(grid, zVal, ...) {
  xyz = cbind(coordinates(grid), grid[[zVal]])
  dimnames(xyz) = list(NULL, c(coordnames(grid), zVal))
  #sort data first by x, then by y, so it fills the matrix correctly...
  xyz = xyz[do.call(order, lapply((2:1), function(i) xyz[, i])), ]
  
  zMatrix = matrix(xyz[,3], grid@grid@cells.dim[1], grid@grid@cells.dim[2])
  SLDF = ContourLines2SLDF(contourLines(x = unique(xyz[,1]), y = unique(xyz[,2]), z = zMatrix, ...), proj4string = grid@proj4string)
  coordnames(SLDF) = coordnames(grid)
  return(SLDF)
}


distAlongSpatialLines = function(x, spLines, maxDist = NA, withAttrs = TRUE) {
  
  if(any(sapply(spLines@lines, function(.l) length(.l@Lines)) > 1)) stop("All Lines in the SpatialLines object 'spLines' must be simplified (e.g., can only contain one sp Line object).  Consider calling explodeSpatialLines() before distanceAlongLines().")
  
  #count the segments in each (simple) lines object
  segCounts = sapply(spLines@lines, function(.l) nrow(.l@Lines[[1]]@coords) - 1)
  #determine the starting value of segment IDs for each Lines object so that all
  #segments have unique ids
  startIDs = c(0, cumsum(segCounts)[-length(segCounts)])
  # get a list of lists of Lines objects, each with one segment of the (simple) Lines
  # objects contained spLines.
  segments = 
    mapply(
      function(.i, .st) explodeLine(spLines@lines[[.i]]@Lines[[1]], .st), 
      1:length(spLines@lines), 
      startIDs,
      SIMPLIFY = F
    )
  
  #determine row IDs for the dataframes, built below.  row IDs are equal to the
  #unique ID's of each segment
  rowIDs = lapply(segments, function(.s) as.character(sapply(.s, function(.s1) .s1@ID)))
  
  #make a SpatialLines object for the segments of each Lines object in spLines  
  segments = lapply(segments, SpatialLines, proj4string = spLines@proj4string)
  # convert SpatialLines objects into SpatialLinesDataFrame objects to store the
  # original Line ID associated with each segment
  segments =
    mapply(
      function(.s, .ID, .rn) 
        SpatialLinesDataFrame(
          .s,
          data = data.frame(
            lineID = rep(.ID, length(.rn)),
            lineSegNum = 1:length(.rn),
            segID = .rn,
            row.names = .rn
          )
        ), 
      segments,
      row.names(spLines),
      rowIDs,
      SIMPLIFY = F
    )
  
  
  
  
  # get the length of each segment
  segLengthList = 
    lapply(
      segments, 
      function(.s) {
        structure(
          sapply(.s@lines, LinesLength),
          names = row.names(.s)
        )
      }
    )
  names(segLengthList) = row.names(spLines)
  # and the cumulative length up to and including each segment, for each line.
  cumLengthList = 
    lapply(
      segLengthList,
      function(.l) {
        structure(
          cumsum(.l),
          names = names(.l)
        )
      }
    )
  names(cumLengthList) = row.names(spLines)
  
  #Because segIDs are unique, we can now dispense with the lists and just make
  #vectors
  segLength = unlist(segLengthList)
  names(segLength) = do.call(c, lapply(segLengthList, names))
  cumLength = unlist(cumLengthList)
  names(cumLength) = names(segLength)
  
  #Make a SpatialLinesDataFrame for all segments across lists.
  segments = do.call(rbind, segments)
  #Snap the points to the closest segment
  x = snapPointsToLines(x, segments, maxDist = maxDist, withAttrs = withAttrs, idField = "segID")
  
  #Pull the segment ID for each point from the results
  segIdx = as.character(x@data$nearest_line_id)
  
  #Get the line IDs associated with each segment
  lineIdx = as.character(segments@data$lineID)
  names(lineIdx) = row.names(segments)
  
  #Now map the line IDs to the points, using segIdx (because there can be more
  #than one point associated with any segment)
  lineIdx = lineIdx[segIdx]
  
  #Get the endpoints of each segment 
  segmentXYs = lapply(segments@lines, function(.l) coordinates(.l@Lines[[1]]))
  names(segmentXYs) = row.names(segments)
  
  #Calculate the fraction of the segment length at the point.
  fracOfSegment = mapply(
    .fractionBetween,
    lapply(1:length(segIdx), function(.i) coordinates(x)[.i,]),
    segmentXYs[segIdx]
  )
  
  #Take the cumsum to the end of the segment; subtract (1 - fracOfSegment)
  distances = cumLength[segIdx] - (segLength[segIdx] * (1 - fracOfSegment))
  
  #adjust the attribute table for the snapped points a bit.
  x@data$lineDistance = distances
  x@data$nearest_line_id = NULL
  x@data$alongLineID = lineIdx
  return(x)
}

#returns a list of Lines objects, each containing a single Line object consisting of one segment
explodeLine = function(x, startID = 0) { #x is a Line objext
  coords = coordinates(x)
  segs = lapply(1:(nrow(coords) - 1), function(.i) Lines(list(Line(coords[.i:(.i+1),])), as.character(.i + startID - 1)))
  return(segs)
}

#returns a list of Lines object, each containing a single Line object of 1 or more segments.
explodeLines = function(x, startID = 0) { #x is a object of type Lines
  return(
    lapply(1:length(x@Lines), function(.i) Lines(list(x@Lines[[.i]]), as.character(.i + startID - 1)))
  )
}

# returns an object of type SpatialLines or SpatialLinesDataFrame (determined by class of x argument)
# with each Line object in x in its own Lines object.
explodeSpatialLines = function(x) {
  lineCounts = sapply(x@lines, function(.x) length(.x@Lines))
  xIDs = data.frame(OrigID = sapply(x@lines, function(.x) return(.x@ID)))
  x@lines = unlist(lapply(x@lines, explodeLines))
  newIDs = 0:(length(x@lines)-1)
  x@lines = mapply(
    function(.x, .ID) {
      .x@ID = as.character(.ID)
      return(.x)
    },
    x@lines,
    newIDs,
    SIMPLIFY = F
  )
  dfIdx = rep(1:nrow(xIDs), lineCounts)
  if(class(x) == "SpatialLinesDataFrame") {
    x@data = structure(
      data.frame(xIDs[dfIdx,], x@data[dfIdx,], row.names = newIDs),
      names = c("OrigID", names(x@data))
    )
  } else {
    x = SpatialLinesDataFrame(x, data.frame(OrigID = xIDs[dfIdx,], row.names = newIDs))
  }
  return(x)
}

.makePattern = function(xLen, patternDistances, offset) {
  patternLength = sum(patternDistances)
  if((patternLength + offset) > xLen) stop("Offset plus sum of patternDistances can't exceed the line length.")
  nRepeat = as.integer((xLen-offset)/patternLength)+1
  patternDistances = c(offset, offset + cumsum(rep(patternDistances, nRepeat)))
  return(patternDistances[patternDistances <= xLen])
}

pointsAlongLines = function(x, distances = fractions * LinesLength(x), fractions = distances / LinesLength(x)) { # x is a Lines object
  if(length(x@Lines) > 1) stop("The Lines object 'x' must be a simple line (e.g., can only have one 'Line' object in the 'Lines' slot).  Consider using 'explode' to convert polylines to simple lines.")
  return(pointsAlongLine(x@Lines[[1]], distances))
}

pointsAlongLine = function(x, distances = fractions * LineLength(x), fractions = distances / LineLength(x)) { # x is a Line object
  if(any(fractions < 0 || fractions > 1.0)) stop("distances argument can not exceed length of line 'x'; fractions argument can not exceed 1.0")
  # get the length and accumulated lengths of each segment
  segLens = segLengths(x)
  cumLengths = cumsum(segLens)
  # find the segement, segment length where the point lies
  targetSeg = sapply(distances, function(.d) min(which(cumLengths >= .d))) #min(which(cumLengths >= distances))
  segLen = segLens[targetSeg]
  endPoints = lapply(targetSeg, function(.t) coordinates(x)[.t:(.t+1),])
  # calculate fractions of the segments to trim
  trimFrac = (cumLengths[targetSeg] - distances)/segLen
  # calculated the point where the line was split
  ptCoord = t(mapply(function(.e, .t) .e[1,] + (.e[2,] - .e[1,]) * (1-.t), endPoints, trimFrac))
  return(ptCoord)
}

repeatAlongLine = function(x, patternDistances, offset = 0, proj4string) {
  return(
    SpatialPointsDataFrame(
      pointsAlongLine(x, .makePattern(LineLength(x), patternDistances, offset)),
      data.frame(distance = patternDistances),
      proj4string = proj4string
    )
  )
}

repeatAlongLines = function(x, patternDistances, offset = 0, proj4string) {
  patternDistances = .makePattern(LinesLength(x), patternDistances, offset)
  return(
    SpatialPointsDataFrame(
      pointsAlongLines(x, patternDistances),
      data.frame(distance = patternDistances),
      proj4string = proj4string
    )
  )
}

segLengths = function(x) {
  coords = coordinates(x)
  nSegs = nrow(coords) - 1
  # pythagorian theorum
  return(sqrt(apply((coords[2:(nSegs +1),] - coords[1:nSegs,])^2, 1, sum)))
}

## This fuction needs to find any vertices where only two lines meet and join 
## the lines into a single line.  When the length of the vector in 
## "connectedLineList" is 1, then only two lines connect at the vertex. Right
## now, the behavior is to connect lines ONLY if ALL of the lines form a single,
## non-networked path.
simplifySpatialLines = function(x, minSpurLength) { #x is Lines object
  x = explodeSpatialLines(x)
  tooShort = sapply(x@lines, LinesLength) < minSpurLength
  
  endPointsList = list(start = lapply(x@lines, function(.x) coordinates(.x)[[1]][1,]))
  endPointsList = c(endPointsList, list(end = lapply(x@lines, function(.x) coordinates(.x)[[1]][nrow(coordinates(.x)[[1]]),])))
  
  # see if first and last point of each line are equal to, repectively, last
  # and first point of every line.  Return index of connected lines.
  connectedLineList = 
    mapply(
      function(i, j)
        lapply(
          endPointsList[[i]], 
          function(ep1)
            which(sapply(endPointsList[[j]], identical, y=ep1))
        ),
      c("start", "end"),
      c("end", "start"),
      SIMPLIFY = F
    )
  # determine which lines are connected to other lines at only one end
  noConnection = sapply(connectedLineList, function(x) lapply(x, length)) == 0
  if(class(noConnection) == "logical") noConnection = matrix(noConnection, nrow = 1)
  spurs = apply(noConnection, 1, function(r) xor(r[1], r[2]))
  # get indexes of lines to keep (those that are not short spurs)
  keepEm = which(!(spurs & tooShort))
  # keep the good lines
  x@lines = x@lines[keepEm]
  # keep track of which are spurs
  spurs = spurs[keepEm]
  # keep only the connection vectors that are relevent
  connectedLineList = lapply(connectedLineList, "[", keepEm)
  # get rid of indices that point to deleted short spurs and update the
  # remaining indices to reflect the new index numbering (resulting from
  # removal of short spurs)
  connectedLineList = lapply(connectedLineList, function(x) lapply(x, function(y) which(keepEm == y[y %in% keepEm])))
  
  if(class(x) == "SpatialLinesDataFrame") x@data = x@data[keepEm,]
  
  # if any line endpoint in connected to more than one other line, the Lines object is still a network.
  isNetwork = any(sapply(connectedLineList, function(x) sapply(x, function(y) length(y)>1)))
  # if there are  two "spurs" left, these are the first and last segments
  # in the polyline.  So if there are 2 spurs and the Lines do not form a
  # network, the Lines can be joined into a single Line.
  
  if(sum(spurs) == 2 & !isNetwork) {
    
    # a little recursive fuction to sort the line based on the connection list   
    sortLines = function(sortedIndexes, count = 1) {
      sortedIndexes[count+1] = connectedLineList$end[[sortedIndexes[count]]]
      count = count + 1
      if(count == length(connectedLineList$end)) {
        return(sortedIndexes)
      } else {
        return(sortLines(sortedIndexes, count))
      }
    }
    
    newOrder = sortLines(which(sapply(connectedLineList$start, length) == 0))
    
    # get order coordinate data.frames
    coordList = unlist(coordinates(x), recursive = F)[newOrder]
    # the last coordinate for all but the last line are duplicates
    dupRows = cumsum(sapply(coordList, nrow))[1:(length(coordList)-1)]
    newCoords = do.call(rbind, coordList)[-dupRows,]
    x@lines = list(Lines(list(Line(newCoords)), "0"))
    if(class(x) == "SpatialLinesDataFrame") x@data = data.frame(ID = 0, row.names = "0")
  }
  return(x)  
}

segmentLine = function(x, segStartDistance, segEndDistance, fractions = F) {
  if (fractions) {
    segStartDistance = segStartDistance * LineLength(x)
    segEndDistance = segEndDistance * LineLength(x)
  }
  if(any(c(segStartDistance, segEndDistance) < 0 || 
         c(segStartDistance, segEndDistance) > LineLength(x))) {
    stop("if fracitons==F, all values in distance argument must be 0 < distances < length(x); if fractions == T, all values must be 0 < fractions < 1.")
  }
  startPoints = pointsAlongLine(x, segStartDistance)
  endPoints = pointsAlongLine(x, segEndDistance)

  vertexDist = c(0, cumsum(segLengths(x)))
#  get the coordinates that go along with the vertexDistances that are between start and end
#  segCoords = lapply(seq(along = startPoints))
}

splitLine = function(x, distances = fractions * LineLength(x), fractions = distances / LineLength(x)) { # x is a Line object
  if(any(fractions<=0.0 || fractions >= 1.0)) stop("all values in distances argument must be 0 < distances < length(x); all values in fractions argument must be 0 < fractions < 1.")
  if(!identical(unique(fractions), fractions)) stop("Distances and fractions can not be duplicated.")
  # get the corrdinates where the line should be split
  splitPts = pointsAlongLine(x, distances)
  # combine with the vertex coordinates
  allCoords = rbind(splitPts, coordinates(x))
  # make a list of the distances for the splits points and vertices
  vertexDist = c(distances, 0, cumsum(segLengths(x)))
  # tag the split points
  isSplitPt = (seq(along = vertexDist) <= length(distances))
  
  # determine the order of the distances
  vertexDistIdx = order(vertexDist)
  
  # reorder everything by distance
  vertexDist = vertexDist[vertexDistIdx]
  allCoords = allCoords[vertexDistIdx,]
  isSplitPt = isSplitPt[vertexDistIdx]
  
  # determine duplicate values
  lVD = length(vertexDist)
  isDuplicate = vertexDist[1:lVD-1] == vertexDist[2:lVD]
  isDuplicate = c(F, isDuplicate) | c(isDuplicate, F)
  
  # Duplicate points are caused with a split point aligns with a line vertex. 
  # Get rid of the duplicated vertex, keeping the one that is tagged as a split
  # point
  keep = isSplitPt | !isDuplicate
  allCoords = allCoords[keep,]
  isSplitPt = isSplitPt[keep]
  
  #make a list of the vertex indices that comprise beginning and end points of
  #the segments.
  endPointIdx = c(1, which(isSplitPt), nrow(allCoords))
  
  LineList = lapply(
    1:(length(endPointIdx) - 1), 
    function (.i) {
      Line(allCoords[endPointIdx[.i]:endPointIdx[.i+1],])
    })
  return(LineList)
}

splitLines2 = function(x, startDistance ) {}

splitLines = function(x, distances = fractions * LinesLength(x), fractions = distances / LinesLength(x), IDs = as.character(seq(0,length(distances)))) { # x is a Lines object
  if(length(x@Lines) > 1) stop("The Lines object 'x' must be a simple line (e.g., can only have one 'Line' object in the 'Lines' slot).  Consider using 'explode' to convert polylines to simple lines.")
  if(length(IDs) != length(distances)+1) stop("The number of ID's must be one greater than the number of distances or fractions provided.")
  if(!identical(unique(IDs), IDs)) stop("IDs must be unique.")
  return(
    structure(
      mapply(
        function(.l, .ID) Lines(.l, .ID),
        splitLine(x@Lines[[1]], distances),
        IDs,
        SIMPLIFY = F
      ),
      names = IDs
    )
  )
}


# # Generate data
# huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
# h <- ggplot(huron, aes(year))
# 
# h + geom_ribbon(aes(ymin=0, ymax=level))
# h + geom_area(aes(y = level))
# 
# # Add aesthetic mappings
# h +
#   geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#   geom_line(aes(y = level))
