# bestOverlaps.R
# 
# Package: t2d
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


bestOverlaps <- function(
  query,
  subject,
  removeSelfMatches           = identical(query, subject)
) {
  allOverlaps <- findOverlaps(query, subject)
  if(removeSelfMatches) {
    if(length(query) != length(subject)) {
      stop("query and subject have different lengths and requested to remove self matches")
    } else {
      allOverlaps <- allOverlaps[queryHits(allOverlaps) != subjectHits(allOverlaps)]
    }
  }
  reciprocalOverlaps <- pmin(
    width(pintersect(query[queryHits(allOverlaps)], subject[subjectHits(allOverlaps)]))/width(query[queryHits(allOverlaps)]),
    width(pintersect(query[queryHits(allOverlaps)], subject[subjectHits(allOverlaps)]))/width(subject[subjectHits(allOverlaps)])
#    width(ranges(allOverlaps, ranges(query), ranges(subject)))/width(query[queryHits(allOverlaps)]),
#    width(ranges(allOverlaps, ranges(query), ranges(subject)))/width(subject[subjectHits(allOverlaps)])
  )
  overlapsForEachQuery <- split(subjectHits(allOverlaps), queryHits(allOverlaps))
  bestMatches <- mapply(
    function(x, y) x[y],
    overlapsForEachQuery,
    tapply(reciprocalOverlaps, queryHits(allOverlaps), which.max)
  )
  matchMatrixMatches <- cbind(
    "query" = as.integer(names(bestMatches)),
    "subject" = bestMatches
  )
  new("Hits",
    queryHits = as.integer(names(bestMatches)),
    subjectHits = bestMatches,
    queryLength = queryLength(allOverlaps),
    subjectLength = subjectLength(allOverlaps)
  )
#  new("RangesMatching",
#    matchMatrix = matchMatrixMatches, 
#    DIM = dim(allOverlaps)
#  )
}

