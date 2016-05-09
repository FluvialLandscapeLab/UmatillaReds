getSegData = function(segBounds, locations, dataAtLocations) {
  segIndexes = lapply(segBounds, function(.b) which((locations>=.b['start']) & (locations<=.b['end'])))
  segLocations = lapply(segIndexes, function(.i) locations[.i])
  segData = lapply(segIndexes, function(.i) dataAtLocations[.i])
  return(list(locations = segLocations, data = segData))
}


calcLongitudinalTrend <- function (location, boundsList, responseVar, pointSpacing, minLength = 0) {
  -unlist(
    mapply(
      localTrend, 
      bounds = boundsList, 
      MoreArgs = list(
        valueVector = responseVar, 
        locationVector = location, 
        pointSpacing = pointSpacing,
        minLength = minLength
      ),
      SIMPLIFY = F
    )
  )
}

steepestSegBounds = function(bounds, locations, dataAtLocations, pointSpacing, posSlope, lookBack) {
  if( (bounds["end"] > (max(locations) + pointSpacing/2)) ||
      (bounds["start"] < (min(locations) - pointSpacing/2)) ) 
  {
    result = bounds
  } else {
    
    segIndexes = which((locations>=bounds['start']) & (locations<=bounds['end']))
    segDataVals = dataAtLocations[segIndexes]
    endDataVal = ifelse(xor(posSlope, lookBack), max(segDataVals), min(segDataVals))
    
    endDataIdx = segIndexes[which(segDataVals == endDataVal)]
    if(lookBack) {
      startDataIdx = max(segIndexes)
      # next line is only necessary in the event that >1 point in the segment has the endDataVal
      endDataIdx = max(endDataIdx)
    } else {
      startDataIdx = min(segIndexes)
      endDataIdx = min(endDataIdx)
    }
    segIndexes = sort(c(startDataIdx, endDataIdx))
    result = structure(locations[segIndexes], names = c("start", "end"))
  }
  return(result)
}

segBounds = function(location, distance, look) {
  relativeDist = list(c(-1,0), c(-0.5, 0.5), c(0,1))
  boundsList = mapply(
    function(.l, .d) structure(.l + .d * relativeDist[[look+2]], names = c("start", "end")),
    location,
    distance,
    SIMPLIFY = F
  )
}

# getStartAndEnd = function(location, distance) {
#   if (distance < 0) {
#     start = location + distance
#     end = location
#   } else {
#     start = location
#     end = location + distance
#   }
#   return(c(start = start, end = end))
# }

localReddCount = function(bounds, redds, samplingPoints) {
  
  if( (samplingPoints[bounds["end"],]$distance > max(redds$lineDistance)) ||
      (samplingPoints[bounds["start"],]$distance < min(redds$lineDistance)) ) 
  {
    result = NA
  } else {
    result = 
      c(
        wide = sum(
          (redds$lineDistance >= samplingPoints[bounds["start"],]$distance) & 
            (redds$lineDistance <= samplingPoints[bounds["trough"],]$distance)
        ),
        narrow = sum(
          (redds$lineDistance >= samplingPoints[bounds["trough"],]$distance) & 
            (redds$lineDistance <= samplingPoints[bounds["end"],]$distance)
        )
      )
  }
  return(result)  
}

localTrend = function(bounds, valueVector, locationVector, pointSpacing = 0, minLength = 0) {
  # location is the location of the downstream point of interest (in same linear units as location vector)
  # distance is the upstream window size
  # valueVector is the response variable
  # locationVector is the line locations corresponding to the values in valueVector
  
  if( (bounds["end"] > (max(locationVector) + pointSpacing/2)) ||
      (bounds["start"] < (min(locationVector) - pointSpacing/2)) ) 
  {
    result = NA
  } else {
    local = which( 
      (locationVector >= (bounds["start"] - pointSpacing/2)) & 
        (locationVector <= (bounds["end"] + pointSpacing/2)) &
        !is.na(valueVector)
    )
    if(length(local) == 0 || length(local)*pointSpacing < minLength) {
      result = NA
    } else {
      model = lm(valueVector[local] ~ locationVector[local])
      result = unname(coefficients(model)[2])
    }
  }
  return(result)
}

calcPeaks = function(x, fun = max, nPreAndPost = 1) {
  xIdx = (1 + nPreAndPost):(length(x) - nPreAndPost)
  local = (-nPreAndPost):nPreAndPost
  is.peak =
    sapply(
      xIdx,
      function(.i) x[.i] == do.call(fun, list(x[.i + local]))
    )
  is.peak = c(rep(NA, nPreAndPost), is.peak, rep(NA, nPreAndPost))
  return(is.peak)
}

pointCounts = function(location, distance, locationVector) {
  
  bounds = getStartAndEnd(location, distance)
  
  pointCount = sum( (locationVector >= bounds["start"]) & (locationVector <= bounds["end"]) )
  return(pointCount)
}

filterDF = function(df, filterList = NULL) {
  #if there is no filter, set the index to true so the whole df is returned
  if(is.null(filterList)) {
    index = T
  } else {
    # create r code text or the "OR" conditions from each sublist
    filterStringList =  lapply(filterList, function(x) paste0("df$", names(x), x, collapse = " | "))
    # parse and evaluate the r code to get boolean vectors
    indexList = lapply(filterStringList, function(filterString) eval(parse(text = filterString)))
    # create a matrix of boolean vectors
    orMatrix = do.call(cbind, indexList)
    # "and" the rows of the matrix by calling all() on each row
    index = apply(orMatrix, 1, all)
  }
  # filter and return the dataframe
  return(df[index,])
}