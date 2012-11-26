# adamsRegionsToExclude.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


adamsRegionsToExclude <- function(
  regionsToDropFilename       = "/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Omni_SVraw/Regions_to_drop_082412.txt"
) {
  adamRegionsToDrop <- scan(regionsToDropFilename, what="")
  adamRegionsToDropDF <- data.frame(
    chrom=sub("^chr([0-9]+):([0-9]+)-([0-9]+)$", "\\1", adamRegionsToDrop),
    start=as.integer(sub("^chr([0-9]+):([0-9]+)-([0-9]+)$", "\\2", adamRegionsToDrop)),
    end=as.integer(sub("^chr([0-9]+):([0-9]+)-([0-9]+)$", "\\3", adamRegionsToDrop)),
    stringsAsFactors=FALSE
  )
  adamRegionsToDropGR <- GRanges(adamRegionsToDropDF[["chrom"]], IRanges(adamRegionsToDropDF[["start"]], adamRegionsToDropDF[["end"]]))
  return(adamRegionsToDropGR)
}
