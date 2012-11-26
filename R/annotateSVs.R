# annotateSVs.R
# 
# Package: t2d
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#library(t2d)
#svsIncNonPASS <- readVcf("/ib/projects/got2d/notBackedUp/tempData/variantCalling/technicalAim918_chr20/sv/all_918bams_chr20run3.sites.vcf", "b37")
#svsGoT2d <- readVcf("/ib/projects/got2d/notBackedUp/tempData/variantCalling/technicalAim918_chr20/sv/918_chr20_pass.genotypes.vcf.gz", "b37")
#svs1kgChr20WithGenotypes <- readVcf("/ib/projects/got2d/notBackedUp/tempData/variantCalling/technicalAim918_chr20/sv/ALL.chr20.phase1_intergrated_calls.20101123.sv.genotypes.vcf.gz", "b37")
#temp <- annotateSVs(svsGoT2d, svs1kgChr20WithGenotypes, annotationStages=c("setEnds", "reciprocalOverlaps"))
#vcf <- annotateSVs(svsGoT2d, svs1kgChr20WithGenotypes)

annotateSVs <- function(
  vcf,
  subjectVcfs                 = NULL,
  annotationStages            = c(
    "setEnds",
    "filterGenotypes",
    "genotypeSummaries",
    "GLsummaries",
    "reciprocalOverlaps",
    "selfReciprocalOverlaps"
#    "filterVariants"
  ),
  endColumnInInfoVcf          = "END",
  endColumnInInfoSubjectVcfs  = paste(names(subjectVcfs), endColumnInInfoVcf, sep="_"),
  vcfGenotypesToAltAlleles    = c("."=-1, "./."=-1, "0/0"=0, "0|0"=0, "0/1"=1, "1/0"=1, "0|1"=1, "1|0"=1, "1/1"=2, "1|1"=2),
  filtersToApply              = c("outOfHWE1e10", "monomorphic", "GSNPAIRS_GSNSAMPLES_RATIO_LT_1.1"),
  verbose                     = TRUE
) {
  require(VariantAnnotation)
  require(snpStats)
  if(!is(vcf, "VCF")) {
    stop("input vcf is not a valid object of class VCF")
  }
  if(!is(subjectVcfs, "list")) {
    stop("subjectVcfs is not a list object")
  }
  
  if("setEnds" %in% annotationStages) {
    if(verbose) cat("annotateSVs: Setting correct SV ends\n------------------------------------\n\n")
    end(ranges(rowData(vcf))) <- values(info(vcf))[[endColumnInInfoVcf]]
    if(!is.null(subjectVcfs)) {
      lapply(
        seq(along=subjectVcfs),
        function(subjectIndex) {
          if(endColumnInInfoSubjectVcfs[subjectIndex] %in% names(values(info(subjectVcfs[[subjectIndex]])))) {
            end(ranges(rowData(subjectVcfs[[subjectIndex]]))) <- values(info(subjectVcfs[[subjectIndex]]))[[endColumnInInfoSubjectVcfs[subjectIndex]]]
          }
        }
      )
    }
    if(verbose) cat("annotateSVs: Creating Lengths\n------------------------------------\n")
    info(vcf) <- BiocGenerics::cbind(values(info(vcf)), DataFrame(Length=width(ranges(rowData(vcf)))-2))
    if(verbose) cat("annotateSVs: Updating vcf header\n\n")
    if("FILTER" %in% names(fixed(exptData(vcf)[["header"]]))) {
      exptData(vcf)[["header"]] <- VCFHeader(
        samples=samples(exptData(vcf)[["header"]]),
        header=DataFrameList(
          META=meta(exptData(vcf)[["header"]]),
          ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
          FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
          FORMAT=geno(exptData(vcf)[["header"]]),
          INFO=BiocGenerics::rbind(
            info(exptData(vcf)[["header"]]),
            DataFrame(Number="1", Type="Integer", Description="Length of deletion", row.names="Length")
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
            DataFrame(Number="1", Type="Integer", Description="Length of deletion", row.names="Length")
          )
        )
      )
    }
  }
  if("filterGenotypes" %in% annotationStages) {
    if(verbose) cat("annotateSVs: Setting genotypes for filtered genotypes to ./.\n------------------------------------\n\n")
    genotypesToFilter <- which(is.na(geno(vcf)[["FT"]]) | geno(vcf)[["FT"]] != "PASS")
#    browser()
#    vcforig <- vcf
    genotypes <- geno(vcf)[["GT"]]
#    class(genotypes)
#    dim(genotypes)
#    table(genotypes, useNA="ifany")
    genotypes[genotypesToFilter] <- "./."
#    which(names(geno(vcf))=="GT")
    genovcf <- geno(vcf)
    genovcf[["GT"]] <- genotypes
#    geno(vcf)[[6]] <- genotypes
#    geno(vcf) <- genovcf
    vcf <- VCF(
      fixed=values(fixed(vcf)),
      info=values(info(vcf)),
      exptData=exptData(vcf),
      rowData=rowData(vcf),
      colData=colData(vcf),
      geno=genovcf
    )
  }
  if("genotypeSummaries" %in% annotationStages) {
    if(verbose) cat("annotateSVs: Creating genotype summaries\n------------------------------------\n")
    vcfGenotypes <- geno(vcf)[["GT"]]
    if(verbose) cat("annotateSVs: Creating SnpMatrix\n")
    genotypesSnpMatrix <- new(
      "SnpMatrix",
      t(
        matrix(
          vcfGenotypesToAltAlleles[vcfGenotypes]+1,
          ncol=ncol(vcfGenotypes),
          dimnames=lapply(dimnames(vcfGenotypes), make.unique)
        )
      )
    )
    if(verbose) cat("annotateSVs: Creating column summary\n")
    genotypesSummary <- col.summary(genotypesSnpMatrix)
    if(verbose) cat("annotateSVs: Analysing ALT alleles\n")
    altAllelesMatrix <- matrix(
      vcfGenotypesToAltAlleles[vcfGenotypes],
      ncol=ncol(vcfGenotypes),
      dimnames=lapply(dimnames(vcfGenotypes), make.unique)
    )
    altAllelesMatrix[altAllelesMatrix==-1] <- NA
    propHetsPerVariant <- apply(
      altAllelesMatrix,
      1,
      function(x) length(which(x==1))/length(!is.na(x))
    )
    if(verbose) cat("annotateSVs: Adding additional columns to genotype summary\n")
    genotypesSummary[["propHetsPerVariant"]] <- propHetsPerVariant
    NAsPerVariant <- apply(
      altAllelesMatrix,
      1,
      function(x) length(which(is.na(x)))
    )
    genotypesSummary[["noCalls"]] <- NAsPerVariant
    genotypesSummary[["allNoCalls"]] <- NAsPerVariant==dim(altAllelesMatrix)[2]
    genotypesSummary[["isMonomorphicRef"]] <- apply(altAllelesMatrix, 1, function(x) length(which(x!=0)) == 0)
    genotypesSummary[["isMonomorphicAlt"]] <- apply(altAllelesMatrix, 1, function(x) length(which(x!=2)) == 0)
    genotypesSummary[["isMonomorphic"]] <- genotypesSummary[["isMonomorphicRef"]] | genotypesSummary[["isMonomorphicAlt"]]
    genotypesSummary[["isOutOfHWE1e10"]] <- genotypesSummary[["z.HWE"]] > -(qnorm(1e-10))
    genotypesSummary[["isOutOfHWE1e10"]][is.na(genotypesSummary[["isOutOfHWE1e10"]])] <- FALSE
    genotypesSummary[["isOutOfHWE0.05"]] <- genotypesSummary[["z.HWE"]] > -(qnorm(0.05))
    genotypesSummary[["isOutOfHWE0.05"]][is.na(genotypesSummary[["isOutOfHWE0.05"]])] <- FALSE
    genotypesSummary[["shouldRemove"]] <- genotypesSummary[["isMonomorphic"]] | genotypesSummary[["isOutOfHWE1e10"]]
#    genotypesSummary[["Length"]] <- width(ranges(rowData(vcf)))-2
    if(verbose) cat("annotateSVs: Adding additional columns to vcf\n")
#    info(vcf) <- BiocGenerics::cbind(values(info(vcf))[, -1], DataFrame(genotypesSummary)) # The -1 avoids duplicate paramRangeID columns - this is a bug in VariantAnnotation that I should report...
    info(vcf) <- BiocGenerics::cbind(values(info(vcf)), DataFrame(genotypesSummary)) # Seems like the bug mentioned above has gone away
    if(verbose) cat("annotateSVs: Updating vcf header\n\n")
    exptData(vcf)[["header"]] <- VCFHeader(
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
        FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=BiocGenerics::rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="1", Type="Integer", Description="snpStats Calls", row.names="Calls"),
          DataFrame(Number="1", Type="Float", Description="snpStats Call.rate", row.names="Call.rate"),
          DataFrame(Number="1", Type="Float", Description="snpStats Certain.calls", row.names="Certain.calls"),
          DataFrame(Number="1", Type="Float", Description="snpStats RAF", row.names="RAF"),
          DataFrame(Number="1", Type="Float", Description="snpStats MAF", row.names="MAF"),
          DataFrame(Number="1", Type="Float", Description="snpStats P.AA", row.names="P.AA"),
          DataFrame(Number="1", Type="Float", Description="snpStats P.AB", row.names="P.AB"),
          DataFrame(Number="1", Type="Float", Description="snpStats P.BB", row.names="P.BB"),
          DataFrame(Number="1", Type="Float", Description="snpStats z.HWE", row.names="z.HWE"),
          DataFrame(Number="1", Type="Float", Description="propHetsPerVariant", row.names="propHetsPerVariant"),
          DataFrame(Number="1", Type="Integer", Description="noCalls", row.names="noCalls"),
          DataFrame(Number="0", Type="Flag", Description="allNoCalls", row.names="allNoCalls"),
          DataFrame(Number="0", Type="Flag", Description="isMonomorphicRef", row.names="isMonomorphicRef"),
          DataFrame(Number="0", Type="Flag", Description="isMonomorphicAlt", row.names="isMonomorphicAlt"),
          DataFrame(Number="0", Type="Flag", Description="isMonomorphic", row.names="isMonomorphic"),
          DataFrame(Number="0", Type="Flag", Description="isOutOfHWE1e10", row.names="isOutOfHWE1e10"),
          DataFrame(Number="0", Type="Flag", Description="isOutOfHWE0.05", row.names="isOutOfHWE0.05"),
          DataFrame(Number="0", Type="Flag", Description="shouldRemove", row.names="shouldRemove")
#          DataFrame(Number="1", Type="Integer", Description="Length of deletion", row.names="Length")
        )
      )
    )
  }
  if("GLsummaries" %in% annotationStages) {
    if(verbose) cat("annotateSVs: Creating GL summaries\n------------------------------------\n")
    GLs <- geno(vcf)[["GL"]]
    if(verbose) cat("annotateSVs: Calculating max GLs\n")
    maxGLs <- apply(GLs, 1, function(x) max(unlist(x)))
    if(verbose) cat("annotateSVs: Adding additional column to vcf\n")
#    info(vcf) <- BiocGenerics::cbind(values(info(vcf))[, -1], DataFrame(maxGL=maxGLs))
    info(vcf) <- BiocGenerics::cbind(values(info(vcf)), DataFrame(maxGL=maxGLs))
    if(verbose) cat("annotateSVs: Updating vcf header\n\n")
    exptData(vcf)[["header"]] <- VCFHeader(
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
        FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=BiocGenerics::rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="1", Type="Float", Description="maximum genotype likelihood", row.names="maxGL")
        )
      )
    )
  }
  if("reciprocalOverlaps" %in% annotationStages) {
    if(is.null(subjectVcfs)) {
      warning("subjectVcfs is NULL, cannot determine reciprocal overlaps")
    } else {
      if(verbose) cat("annotateSVs: Determining reciprocal overlaps and merging\n---------------------------------------------------------\n\n")
      lapply(
        seq(along=subjectVcfs),
        function(subjectIndex) {
          vcf <<- attachBestOverlaps(vcf, subjectVcfs[[subjectIndex]], fieldPrefixForSubjectVcf=paste(names(subjectVcfs)[subjectIndex], "_", sep=""))
        }
      )

    }
  }
  if("selfReciprocalOverlaps" %in% annotationStages) {
    vcf <<- attachBestOverlaps(vcf, vcf, infoFieldsToIncludeFromSubj=NULL, fieldPrefixForSubjectVcf="self_")
  }
  if("filterVariants" %in% annotationStages) {
    filters <- values(fixed(vcf))[["FILTER"]]
#    filters <- character(dim(vcf)[1])
#    filters <- NA
#    filters <- as.character(rep(NA, dim(vcf)[1]))
    if("outOfHWE1e10" %in% filtersToApply) {
      filters[which(values(info(vcf))[["z.HWE"]] > -(qnorm(1e-10)))] <- paste(filters[which(values(info(vcf))[["z.HWE"]] > -(qnorm(1e-10)))], "outOfHWE1e10", sep=";")
    }
    if("monomorphic" %in% filtersToApply) {
      filters[which(values(info(vcf))[["isMonomorphic"]])] <- paste(filters[which(values(info(vcf))[["isMonomorphic"]])], "monomorphic", sep=";")
    }
    if("GSNPAIRS_GSNSAMPLES_RATIO_LT_1.1" %in% filtersToApply) {
      filters[which(values(info(vcf))[["GSNPAIRS"]] / values(info(vcf))[["GSNSAMPLES"]] <= 1.1)] <- paste(filters[which(values(info(vcf))[["GSNPAIRS"]] / values(info(vcf))[["GSNSAMPLES"]] <= 1.1)], "GSNPAIRS_GSNSAMPLES_RATIO_LT_1.1", sep=";")
    }
    filters[filters != "PASS"] <- sub("^PASS;", "", filters[filters != "PASS"])
    fixedVcf <- fixed(vcf)[, -1]
    values(fixedVcf)[["FILTER"]] <- filters
    fixed(vcf) <- values(fixedVcf)
    headerFilters <- BiocGenerics::rbind(
      fixed(exptData(vcf)[["header"]])[["FILTER"]],
      DataFrame(Description="Hardy-Weinberg equilbrium p-value (calculated using hard genotype calls) < 1e-10", row.names="outOfHWE1e10"),
      DataFrame(Description="Variant has monomorphic hard genotype calls in given samples", row.names="monomorphic"),
      DataFrame(Description="GSNPAIRS/GSNSAMPLES <= 1.1", row.names="GSNPAIRS_GSNSAMPLES_RATIO_LT_1.1")
    )
    exptData(vcf)[["header"]] <- VCFHeader(
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
        FILTER=headerFilters,
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=info(exptData(vcf)[["header"]])
      )
    )
#    fixed(exptData(vcf)[["header"]])[["FILTER"]] <- headerFilters
  }
  if(verbose) cat("annotateSVs: Done\n\n")
  return(vcf)
}
