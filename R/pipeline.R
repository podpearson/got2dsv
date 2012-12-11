# pipeline.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# X chromosome run
# vcfX <- pipeline(chromosomes="X", outputVcfRdaFilename="vcfX.rda", outputVcfListRdaFilename="vcfListX.rda")

pipeline <- function(
  outputDir                   = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation",
  outputVcfListRdaFilename    = "vcfList.rda",
  outputVcfRdaFilename        = "vcf.rda",
  shouldReload                = !file.exists(file.path(outputDir, outputVcfRdaFilename)),
  chromosomes                 = c(1:22),
  shouldCreateFilesForHyun    = !file.exists(file.path(outputDir, sprintf("t2dgo_chr%s_stg1_merged.genotypes.fixed.annotatedForHyun.vcf", chromosomes[length(chromosomes)])))
) {
  if(shouldReload) {
#  vcfOmni <- readOmniVcfToR("/ddn/projects11/got2d/rpearson/SVGtypes/GoT2D.Omni.Deletions.092512.vcf", outputDir=outputDir)
#  vcfOmniRaw <- readOmniVcfToR("/ddn/projects11/got2d/rpearson/SVGtypes/GoT2D.Omni.Raw.Deletions.101512_v1.vcf", outputDir=outputDir)
    vcf1kg <- read1kgVcfToR(
      "/ddn/projects11/got2d/GoT2DSVs/1000G_SVgenotypes/ALL.genome.phase1_release_v3.20101123.svs.sites.vcf",
      outputDir=outputDir,
      shouldReload=shouldReload
    )
    vcf1kgUsedInDiscovery <- read1kgDiscoveryVcfToR(
      "/ddn/projects11/got2d/GoT2DSVs/1000G_SVgenotypes/ALL.wgs.merged_5_del_call_sets_bps.20101123.sv_dels.low_coverage.sites.vcf",
      outputDir=outputDir,
      shouldReload=shouldReload
    )
    regionsToExclude <- adamsRegionsToExclude("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Omni_SVraw/Regions_to_drop_082412.txt")
    bobHighVPS <- read.delim("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/bobHighLowVariantCountsEmail20121016/high_vps_samples.dat", as.is=TRUE)[["SAMPLE"]]
  
    vcfList <- sapply(                                                          # this reads is data from relevant vcf files into a VCF object
      chromosomes,                                                              # loops x over the values 1 to 22, runs loadAndAnnotateLowpassSVs for each value of x, and creates a list object containing the 22 results. Alternatively, when run with chromosomes="X" will create a list object with only the result for X chromosome.
      function(x) {
        loadAndAnnotateLowpassSVs(
          chromosome                  = x,
          vcf1kg                      = vcf1kg,
          vcf1kgUsedInDiscovery       = vcf1kgUsedInDiscovery,
          regionsToExclude            = regionsToExclude,
          bobHighVPS                  = bobHighVPS,
          shouldReload                = TRUE
  #        shouldReload                = shouldReload
        )
      }
    )
    save(vcfList, file=file.path(outputDir, outputVcfListRdaFilename))                     # save as an R object for quicker reading in later
    if(length(vcfList) > 1) {
      vcf <- combineVcfListIntoVcf(vcfList)                                       # create a single VCF object containing all chromosomes for use in subsequent analysis. Note this is very memory intensive (peaks at > 16GB) for reasons I don't understand
    } else {
      vcf <- vcfList[[1]]
    }
    gc()                                                                        # this is garabge collection to help keep RAM usage down
    save(vcf, file=file.path(outputDir, outputVcfRdaFilename))                             # save as an R object for quicker reading in later
  } else {
    load(file.path(outputDir, outputVcfRdaFilename))
  }
  if(shouldCreateFilesForHyun) {
    vcfListImportantColumns <- lapply(vcfList, extraImportantInfo)              # return a list object containing VCF objects which have only the information required by Hyun in them. Essentially calls extraImportantInfo function on each element of vcfList, and returns a new list of these
    lapply(                                                                     # writes out to actual vcf files. These are the files I emailed Hyun about 17/11/2012. The index=TRUE ensures output files and bgzipped and tabix indexed
      chromosomes,
      function(chromosome) {
        cat(".")
        writeVcf(
          vcfListImportantColumns[[chromosome]],
          filename = file.path(outputDir, sprintf("t2dgo_chr%s_stg1_merged.genotypes.fixed.annotatedForHyun.vcf", chromosome)),
          index = TRUE
        )
      }
    )
  }
  
##  Some debugging stuff - remove later
#    thunderVcf <- readVcf(
#      "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation/got2d.2874.chr1.sv.thunder.vcf",
#      genome="dummy",
#      param=ScanVcfParam(fixed="FILTER", info=c("END"), geno=NA)
##      param=ScanVcfParam(fixed="FILTER", info=c("AC", "AN", "END", "GCLENGTH", "GLALTFREQ", "SVLEN", "BAVGPOST", "BRSQ", "LDAF", "AVGPOST", "RSQ"), geno=c("GD"))
#    )
#    end(ranges(thunderVcf)) <- values(info(thunderVcf))[["END"]]

#  return(vcfList)
  return(vcf)
}

#lowpassInfo <- do.call(
#  c,
#  lapply(
#    vcfListImportantColumns,
#    function(vcf) {
#      cat(".")
#      temp <- info(vcf)
#      if("set" %in% names(values(temp))) {
#        values(temp)[["set"]] <- NULL
#      }
#      temp
#    }
#  )
#)
#lowpassFixed <- do.call(
#  c,
#  lapply(
#    vcfListImportantColumns,
#    function(vcf) {
#      cat(".")
#      fixed(vcf)
#    }
#  )
#)
#
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(GLALTFREQ==0, sapply(CalledBy, function(x) "GenomeStrip" %in% x), useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(GLALTFREQ==0, paste(ROgt0.8withBatch1, ROgt0.8withBatch2, ROgt0.8with1000G), useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(unlist(SVLEN)==0, sapply(CalledBy, function(x) "GenomeStrip" %in% x), useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(sapply(CalledBy, function(x) paste(x, collapse=":")), unlist(SVLEN)==0, useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(sapply(CalledBy, function(x) paste(x, collapse=":")), is.na(unlist(SVLEN)), useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(is.na(unlist(SVLEN)), batch, useNA="ifany"))
#with(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]), table(unlist(SVLEN)==0, paste(ROgt0.8withBatch1, ROgt0.8withBatch2, ROgt0.8with1000G), useNA="ifany"))
#table(sapply(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"])[["SVLEN"]], length), useNA="ifany")
#
#SVLENs <- unlist(values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"])[["SVLEN"]])
#ENDs <- values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"])[["END"]]
#POSs <- start(ranges(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"]))
#LENGTHS <- POSs-ENDs
#GSELENGTHs <- values(lowpassInfo[values(lowpassFixed)[["FILTER"]]=="PASS"])[["GSELENGTH"]]
#table(GSELENGTHs-(-LENGTHS), useNA="ifany")
#
#chr1info <- values(info(vcfList[[1]]))
#chr1filter <- filt(vcfList[[1]])
#chr1infoPASS <- chr1info[chr1filter=="PASS", ]


