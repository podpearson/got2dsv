# attachBestOverlaps.R
# 
# Package: t2d
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


attachBestOverlaps <- function(
  vcf,
  subjectVcf,
  infoFieldsToIncludeFromSubj = names(values(info(subjectVcf))),
#  infoFieldsToIncludeFromSubj = names(values(info(subjectVcf)))[-1],
#  fieldPrefixForSubjectVcf    = "bestOverlap_",
  fieldPrefixForSubjectVcf    = "subject_",
  removeSelfMatches           = identical(vcf, subjectVcf),
  verbose                     = TRUE
) {
  if(verbose) cat("attachBestOverlaps:\n-------------------\nQuery:\n------\n\n")
  if(verbose) print(vcf)
  if(verbose) cat("\nSubject:\n--------\n\n")
  if(verbose) print(subjectVcf)
  if(verbose) cat("\n")
  query <- rowData(vcf)
  subject <- rowData(subjectVcf)
  allOverlaps <- findOverlaps(query, subject)
  if(removeSelfMatches) {
    if(dim(vcf)[1] != dim(subjectVcf)[1]) {
      stop("vcf and subjectVcf have different numbers of variants and requested to remove self matches")
    } else {
      allOverlaps <- allOverlaps[queryHits(allOverlaps) != subjectHits(allOverlaps)]
    }
  }
  bestReciprocalOverlaps <- bestOverlaps(query, subject, removeSelfMatches)
  numbersOfOverlaps <- table(queryHits(allOverlaps))
  bestOverlapsDF <- DataFrame(numberOfOverlaps=integer(length(query)))
  bestOverlapsDF[["numberOfOverlaps"]] <- NULL
#  bestOverlapsDF <- DataFrame()
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")]] <- integer(length(query))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")]][as.integer(names(numbersOfOverlaps))] <- numbersOfOverlaps
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Name", sep="")]] <- character(length(query))
  if(is.null(names(ranges(subject[subjectHits(bestReciprocalOverlaps)])))) {
    bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Name", sep="")]][queryHits(bestReciprocalOverlaps)] <- NA
  } else {
    bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Name", sep="")]][queryHits(bestReciprocalOverlaps)] <- names(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
  }
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Start", sep="")]] <- integer(length(query))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Start", sep="")]][queryHits(bestReciprocalOverlaps)] <- start(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "End", sep="")]] <- integer(length(query))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "End", sep="")]][queryHits(bestReciprocalOverlaps)] <- end(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Filter", sep="")]] <- character(length(query))
  if("FILTER" %in% names(values(fixed(subjectVcf)))) {
    bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Filter", sep="")]][queryHits(bestReciprocalOverlaps)] <- values(fixed(subjectVcf))[subjectHits(bestReciprocalOverlaps), "FILTER"]
  } else {
    bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "Filter", sep="")]][queryHits(bestReciprocalOverlaps)] <- NA
  }
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep="")]] <- numeric(length(query))
  bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep="")]][queryHits(bestReciprocalOverlaps)]  <- pmin(
    width(pintersect(query[queryHits(bestReciprocalOverlaps)], subject[subjectHits(bestReciprocalOverlaps)]))/width(query[queryHits(bestReciprocalOverlaps)]),
    width(pintersect(query[queryHits(bestReciprocalOverlaps)], subject[subjectHits(bestReciprocalOverlaps)]))/width(subject[subjectHits(bestReciprocalOverlaps)])
  )

#  bestOverlapsDF[["bestOverlapFILTER"]][queryHits(bestReciprocalOverlaps)] <- values(subject[subjectHits(bestReciprocalOverlaps)])[["FILTER"]]

  if(!is.null(infoFieldsToIncludeFromSubj)) {
    tempSubjectInfo <- DataFrame(
      do.call(
        BiocGenerics::rbind,
        lapply(
          seq(length=queryLength(bestReciprocalOverlaps)),
          function(x) {
            values(info(subjectVcf))[1, infoFieldsToIncludeFromSubj]
          }
        )
      )
    )
    names(tempSubjectInfo) <- infoFieldsToIncludeFromSubj
#    dimnames(tempSubjectInfo) <- list(dimnames(tempSubjectInfo)[[1]], infoFieldsToIncludeFromSubj)
    tempSubjectInfo[queryHits(bestReciprocalOverlaps), infoFieldsToIncludeFromSubj] <- values(info(subjectVcf))[subjectHits(bestReciprocalOverlaps), infoFieldsToIncludeFromSubj]
    rowsNotInSubjectVcf <- setdiff(seq(length=queryLength(bestReciprocalOverlaps)), queryHits(bestReciprocalOverlaps))
    tempSubjectInfo[rowsNotInSubjectVcf, infoFieldsToIncludeFromSubj] <- NA
    bestOverlapsDF <- BiocGenerics::cbind(
      bestOverlapsDF,
      tempSubjectInfo
    )
  }
  
#  lapply(
#    infoFieldsToIncludeFromSubj,
#    function(fieldName) {
#      tempSubjectInfo[[fieldName]] <<- NA
#    }
#  )
#  
#  bestOverlapsDF[queryHits(bestReciprocalOverlaps), ] <- BiocGenerics::cbind(
#    bestOverlapsDF[queryHits(bestReciprocalOverlaps), ],
#    values(info(subjectVcf))[subjectHits(bestReciprocalOverlaps), infoFieldsToIncludeFromSubj]
#  )
#  temp <- integer(queryLength(bestReciprocalOverlaps))
#  temp <- rep(NA, queryLength(bestReciprocalOverlaps))
#  temp[queryHits(bestReciprocalOverlaps)] <- subjectHits(bestReciprocalOverlaps)
#  bestOverlapsDF <- BiocGenerics::cbind(
#    bestOverlapsDF,
#    values(info(subjectVcf))[temp, infoFieldsToIncludeFromSubj]
#  )
#  
#  bestOverlapsDF[queryHits(bestReciprocalOverlaps), infoFieldsToIncludeFromSubj] <- values(info(subjectVcf))[subjectHits(bestReciprocalOverlaps), infoFieldsToIncludeFromSubj]
#  
#  
#  lapply(
#    infoFieldsToIncludeFromSubj,
#    function(infoField) {
#      if(verbose) {cat(infoField, "")}
#      infoFieldClass <- class(values(info(subjectVcf))[[infoField]])
#      if(infoFieldClass == "array") {
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]] <- array(dim=c(length(query), dim(info(subjectVcf)[[infoField]])[-1]))
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]][as.integer(names(numbersOfOverlaps)),,] <- values(info(subjectVcf))[[infoField]][subjectHits(bestReciprocalOverlaps),,]
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- array(dim=c(length(query), dim(info(subjectVcf)[[infoField]])[-1]))
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps)),,] <<- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps),,]
###        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps)),,,drop=FALSE] <- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps),,,drop=FALSE]
#      } else if(infoFieldClass == "CompressedIntegerList") {
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]] <- integer(length(query))
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]][as.integer(names(numbersOfOverlaps)),,] <- unlist(values(info(subjectVcf))[[infoField]][subjectHits(bestReciprocalOverlaps)])
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- integer(length(query))
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- unlist(info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)])
#      } else if(infoFieldClass == "CompressedCharacterList") {
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]] <- character(length(query))
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]][as.integer(names(numbersOfOverlaps)),,] <- unlist(values(info(subjectVcf))[[infoField]][subjectHits(bestReciprocalOverlaps)])
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- character(length(query))
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- unlist(info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)])
#      } else {
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]] <- do.call(BiocGenerics::sapply(values(info(subjectVcf)), class)[infoField], list(length(query)))
#        bestOverlapsDF[[paste(fieldPrefixForSubjectVcf, infoField, sep="")]][as.integer(names(numbersOfOverlaps)),,] <- values(info(subjectVcf))[[infoField]][subjectHits(bestReciprocalOverlaps)]
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- do.call(BiocGenerics::sapply(info(subjectVcf), class)[infoField], list(length(query)))
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)]
#      }
#    }
#  )
  
  
#  info(vcf) <- BiocGenerics::cbind(values(info(vcf))[, -1], bestOverlapsDF) # The -1 avoids duplicate paramRangeID columns - this is a bug in VariantAnnotation that I should report...  
  info(vcf) <- BiocGenerics::cbind(values(info(vcf)), bestOverlapsDF) # Seems like the bug mentioned above has gone away
  
  
  if(verbose) cat("attachBestOverlaps: Updating vcf header\n\n")
  if(is.null(fixed(exptData(vcf)[["header"]])[["FILTER"]])) {
    if(!is.null(infoFieldsToIncludeFromSubj)) {
      exptData(vcf)[["header"]] <- VCFHeader(
        samples=samples(exptData(vcf)[["header"]]),
        header=DataFrameList(
          META=meta(exptData(vcf)[["header"]]),
          ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
          FORMAT=geno(exptData(vcf)[["header"]]),
          INFO=BiocGenerics::rbind(
            info(exptData(vcf)[["header"]]),
            DataFrame(Number="1", Type="Integer", Description="number of overlaps in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Name of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Name", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Start coordinate of best overlap in subject vcf (0 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "Start", sep="")),
            DataFrame(Number="1", Type="Integer", Description="End coordinate of best overlap in subject vcf (-1 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "End", sep="")),
            DataFrame(Number="1", Type="String", Description="FILTER status of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Filter", sep="")),
            DataFrame(Number="1", Type="Float", Description="Reciprocal overlap with best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep="")),
            info(exptData(subjectVcf)[["header"]])
          )
        )
      )
    } else {
      exptData(vcf)[["header"]] <- VCFHeader(
        samples=samples(exptData(vcf)[["header"]]),
        header=DataFrameList(
          META=meta(exptData(vcf)[["header"]]),
          ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
          FORMAT=geno(exptData(vcf)[["header"]]),
          INFO=BiocGenerics::rbind(
            info(exptData(vcf)[["header"]]),
            DataFrame(Number="1", Type="Integer", Description="number of overlaps in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Name of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Name", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Start coordinate of best overlap in subject vcf (0 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "Start", sep="")),
            DataFrame(Number="1", Type="Integer", Description="End coordinate of best overlap in subject vcf (-1 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "End", sep="")),
            DataFrame(Number="1", Type="String", Description="FILTER status of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Filter", sep="")),
            DataFrame(Number="1", Type="Float", Description="Reciprocal overlap with best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep=""))
          )
        )
      )
    }
  } else {
    if(!is.null(infoFieldsToIncludeFromSubj)) {
      exptData(vcf)[["header"]] <- VCFHeader(
        samples=samples(exptData(vcf)[["header"]]),
        header=DataFrameList(
          META=meta(exptData(vcf)[["header"]]),
          ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
          FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
          FORMAT=geno(exptData(vcf)[["header"]]),
          INFO=BiocGenerics::rbind(
            info(exptData(vcf)[["header"]]),
            DataFrame(Number="1", Type="Integer", Description="number of overlaps in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Name of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Name", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Start coordinate of best overlap in subject vcf (0 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "Start", sep="")),
            DataFrame(Number="1", Type="Integer", Description="End coordinate of best overlap in subject vcf (-1 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "End", sep="")),
            DataFrame(Number="1", Type="String", Description="FILTER status of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Filter", sep="")),
            DataFrame(Number="1", Type="Float", Description="Reciprocal overlap with best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep="")),
            info(exptData(subjectVcf)[["header"]])
          )
        )
      )
    } else {
      exptData(vcf)[["header"]] <- VCFHeader(
        samples=samples(exptData(vcf)[["header"]]),
        header=DataFrameList(
          META=meta(exptData(vcf)[["header"]]),
          ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
          FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
          FORMAT=geno(exptData(vcf)[["header"]]),
          INFO=BiocGenerics::rbind(
            info(exptData(vcf)[["header"]]),
            DataFrame(Number="1", Type="Integer", Description="number of overlaps in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "NumberOfOverlaps", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Name of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Name", sep="")),
            DataFrame(Number="1", Type="Integer", Description="Start coordinate of best overlap in subject vcf (0 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "Start", sep="")),
            DataFrame(Number="1", Type="Integer", Description="End coordinate of best overlap in subject vcf (-1 means no overlap)", row.names=paste(fieldPrefixForSubjectVcf, "End", sep="")),
            DataFrame(Number="1", Type="String", Description="FILTER status of best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "Filter", sep="")),
            DataFrame(Number="1", Type="Float", Description="Reciprocal overlap with best overlap in subject vcf", row.names=paste(fieldPrefixForSubjectVcf, "ReciprocalOverlap", sep=""))
          )
        )
      )
    }
    
  }

#  info(vcf)[["numberOfOverlaps"]] <- integer(length(query))
#  info(vcf)[["numberOfOverlaps"]][as.integer(names(numbersOfOverlaps))] <- numbersOfOverlaps
#  exptData(vcf)[["HEADER"]][["INFO"]] <- BiocGenerics::rbind(
#    exptData(vcf)[["HEADER"]][["INFO"]],
#    DataFrame(Number="1", Type="Integer", Description="number of overlaps in subject vcf", row.names="numberOfOverlaps")
#  )
##  values(query)[["numberOfOverlaps"]] <- integer(length(query))
##  values(query)[
##    as.integer(names(numbersOfOverlaps)),
##    "numberOfOverlaps"
##  ] <- numbersOfOverlaps
#  info(vcf)[["bestOverlapName"]] <- character(length(query))
#  info(vcf)[["bestOverlapName"]][queryHits(bestReciprocalOverlaps)] <- names(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
#  info(vcf)[["bestOverlapStart"]] <- integer(length(query))
#  info(vcf)[["bestOverlapStart"]][queryHits(bestReciprocalOverlaps)] <- start(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
#  info(vcf)[["bestOverlapEnd"]] <- integer(length(query))
#  info(vcf)[["bestOverlapEnd"]][queryHits(bestReciprocalOverlaps)] <- end(ranges(subject[subjectHits(bestReciprocalOverlaps)]))
#  info(vcf)[["bestOverlapFILTER"]] <- character(length(query))
#  info(vcf)[["bestOverlapFILTER"]][queryHits(bestReciprocalOverlaps)] <- values(subject[subjectHits(bestReciprocalOverlaps)])[["FILTER"]]
##  info(vcf)[["bestOverlapRanges"]] <- new(
##    "IRanges",
##    start=as.integer(rep(0, length(query))),
##    width=as.integer(rep(0, length(query)))
##  )
##  info(vcf)[["bestOverlapRanges"]][queryHits(bestReciprocalOverlaps)] <- ranges(subject[subjectHits(bestReciprocalOverlaps)])
#  exptData(vcf)[["HEADER"]][["INFO"]] <- BiocGenerics::rbind(
#    exptData(vcf)[["HEADER"]][["INFO"]],
#    DataFrame(Number="1", Type="Integer", Description="Start coordinate of best overlap in subject vcf (0 means no overlap)", row.names="bestOverlapStart"),
#    DataFrame(Number="1", Type="Integer", Description="End coordinate of best overlap in subject vcf (-1 means no overlap)", row.names="bestOverlapEnd"),
#    DataFrame(Number="1", Type="String", Description="FILTER status of best overlap in subject vcf", row.names="bestOverlapFILTER")
#  )
#  exptData(vcf)[["HEADER"]][["INFO"]] <- BiocGenerics::rbind(
#    exptData(vcf)[["HEADER"]][["INFO"]],
#    DataFrame(Number="2", Type="Integer", Description="range of best overlap in subject vcf", row.names="bestOverlapRanges")
#  )
#  values(query)[["bestOverlapRanges"]] <- new(
#    "IRanges",
#    start=as.integer(rep(0, length(query))),
#    width=as.integer(rep(0, length(query)))
#  )
#  values(query)[
#    queryHits(bestReciprocalOverlaps),
#    "bestOverlapRanges"
#  ] <- ranges(subject[subjectHits(bestReciprocalOverlaps)])
#  if(verbose) cat("attachBestOverlaps: Merging results into vcf\n")
#  lapply(
#    infoFieldsToIncludeFromSubj,
#    function(infoField) {
#      if(verbose) {cat(infoField, "")}
#      infoFieldClass <- class(info(subjectVcf)[[infoField]])
#      if(infoFieldClass == "array") {
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- array(dim=c(length(query), dim(info(subjectVcf)[[infoField]])[-1]))
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps)),,] <<- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps),,]
##        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps)),,,drop=FALSE] <- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps),,,drop=FALSE]
#      } else if(infoFieldClass == "CompressedIntegerList") {
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- integer(length(query))
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- unlist(info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)])
#      } else if(infoFieldClass == "CompressedCharacterList") {
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- character(length(query))
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- unlist(info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)])
#      } else {
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]] <<- do.call(BiocGenerics::sapply(info(subjectVcf), class)[infoField], list(length(query)))
#        info(vcf)[[paste("bestOverlap", infoField, sep="_")]][as.integer(names(numbersOfOverlaps))] <<- info(subjectVcf)[[infoField]][subjectHits(bestReciprocalOverlaps)]
#      }
##      values(query)[
##        queryHits(bestReciprocalOverlaps),
##        paste("bestOverlap", columnName, sep="_")
##      ] <<- values(subject) [
##                  subjectHits(bestReciprocalOverlaps),
##                  columnName
##              ]
#    }
#  )
#  if(verbose) cat("\nattachBestOverlaps: Updating vcf header\n")
#  subjectHeaderInfo <- exptData(subjectVcf)[["HEADER"]][["INFO"]]
#  subjectHeaderInfo <- subjectHeaderInfo[infoFieldsToIncludeFromSubj, ]
#  row.names(subjectHeaderInfo) <- paste("bestOverlap", row.names(subjectHeaderInfo), sep="_")
#  exptData(vcf)[["HEADER"]][["INFO"]] <- BiocGenerics::rbind(
#    exptData(vcf)[["HEADER"]][["INFO"]],
#    subjectHeaderInfo
#  )
  
#  lapply(
#    names(values(subject)),
#    function(columnName) {
#      if(sapply(values(subject), class)[columnName] %in% "DNAStringSet")
#      values(query)[[paste("bestOverlap", columnName, sep="_")]] <<- do.call(sapply(values(subject), class)[columnName], list(length(query)))
#      values(query)[
#        queryHits(bestReciprocalOverlaps),
#        paste("bestOverlap", columnName, sep="_")
#      ] <<- values(subject) [
#                  subjectHits(bestReciprocalOverlaps),
#                  columnName
#              ]
#    }
#  )
#  if(verbose) cat("attachBestOverlaps: Adding reciprocal overlaps\n")
#  reciprocalOverlaps <- pmin(
#    width(pintersect(query[queryHits(bestReciprocalOverlaps)], subject[subjectHits(bestReciprocalOverlaps)]))/width(query[queryHits(bestReciprocalOverlaps)]),
#    width(pintersect(query[queryHits(bestReciprocalOverlaps)], subject[subjectHits(bestReciprocalOverlaps)]))/width(subject[subjectHits(bestReciprocalOverlaps)])
##    width(ranges(bestReciprocalOverlaps, ranges(query), ranges(subject)))/width(query[queryHits(bestReciprocalOverlaps)]),
##    width(ranges(bestReciprocalOverlaps, ranges(query), ranges(subject)))/width(subject[subjectHits(bestReciprocalOverlaps)])
#  )
#  info(vcf)[["reciprocalOverlap"]] <- numeric(length(query))
#  info(vcf)[["reciprocalOverlap"]][queryHits(bestReciprocalOverlaps)] <- reciprocalOverlaps
#  exptData(vcf)[["HEADER"]][["INFO"]] <- BiocGenerics::rbind(
#    exptData(vcf)[["HEADER"]][["INFO"]],
#    DataFrame(Number="1", Type="Float", Description="reciprocal overlap of best overlap in subject vcf", row.names="reciprocalOverlap")
#  )
#  values(query)[["reciprocalOverlap"]] <- numeric(length(query))
#  values(query)[
#    queryHits(bestReciprocalOverlaps),
#    "reciprocalOverlap"
#  ] <- reciprocalOverlaps
  
  if(verbose) cat("attachBestOverlaps: Done\n\n")
  
  return(vcf)
}

