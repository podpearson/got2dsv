# loadAndAnnotateLowpassSVs.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndAnnotateLowpassSVs <- function(
  chromosome                  = 20,
  got2dVcfFmtBeforeFixing     = "/ddn/projects11/got2d/rpearson/SVGtypes/t2dgo_chr%s_stg1_merged.genotypes.vcf.gz",  # Note the %s as this is used in sprintf to combine with chromsome
  got2dVcfFmt                 = "/ddn/projects11/got2d/rpearson/SVGtypes/SVG9oct/t2dgo_chr%s_stg1_merged.genotypes.fixed.vcf",  # Note the %s as this is used in sprintf to combine with chromsome
  thunderVcfFmt               = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation/got2d.2874.chr%s.sv.thunder.vcf",  # Note the %s as this is used in sprintf to combine with chromsome
#  thunderVcfFmt               = "/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Dwnldd/INTGfinal/svonly/got2d.2874.chr%s.sv.thunder.vcf.gz",
  batch1Fmt                   = "/ddn/projects11/got2d/GoT2DSVs/discVCFmerge/PASS/1510_chr%s_stg1_6filtered_PASS.discovery.vcf",  # Note the %s as this is used in sprintf to combine with chromsome
  batch2Fmt                   = "/ddn/projects11/got2d/GoT2DSVs/discVCFmerge/PASS/1292_chr%s_stg1_6filtered_PASS.discovery.vcf",  # Note the %s as this is used in sprintf to combine with chromsome
  vcfOmni                     = NULL,                                           # not used here
  vcf1kg                      = NULL,                                           # should be passed in by the calling function (pipeline)
  vcf1kgUsedInDiscovery       = NULL,                                           # should be passed in by the calling function (pipeline)
  regionsToExclude            = adamsRegionsToExclude("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Omni_SVraw/Regions_to_drop_082412.txt"),    # here the default values is a call to another function
  bobHighVPS                  = read.delim("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/bobHighLowVariantCountsEmail20121016/high_vps_samples.dat", as.is=TRUE)[["SAMPLE"]],        # here the default values is a call to another function
  outputDir                   = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation",
  outputFilenamePreAnnotation = file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "rda", sep=".")),
  outputFilenamePostAnnotation= file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "annotateSVs.rda", sep=".")),
  outputVcfPostAnnotation     = file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "annotateSVs.vcf", sep=".")),
  outputFilenameThunder       = file.path(outputDir, paste(sprintf(basename(thunderVcfFmt), chromosome), "rda", sep=".")),
  outputFilenameBatch1        = file.path(outputDir, paste(sprintf(basename(batch1Fmt), chromosome), "rda", sep=".")),
  outputFilenameBatch2        = file.path(outputDir, paste(sprintf(basename(batch2Fmt), chromosome), "rda", sep=".")),
  shouldReload                = !file.exists(outputFilenamePostAnnotation),
  shouldAnnotateWithBatchNum  = file.exists(sprintf(batch1Fmt, chromosome)) & file.exists(sprintf(batch2Fmt, chromosome))
) {
  if(shouldReload) {                                                            # have set this up to ensure saved R object files are used if previously created, saving a lot of time when rerunning/debugging
#    the following commands were used to clean up previous versions of vcfs. Not needed here but kept in for reference
#    cat("Adding extra info lines to vcf\n")
#    bashCommand <- paste(
#      "(cat /ddn/projects11/got2d/GoT2DSVs/SVGtypes/Dwnldd/t2dgo_chr", chromosome, "_stg1_merged.genotypes.vcf | head -1000 | grep ^##; ",
#      "echo '##INFO=<ID=GCFRACTION,Number=1,Type=Float,Description=\"GC content fraction\">'; ",
#      "echo '##INFO=<ID=CalledBy,Number=.,Type=String,Description=\"SV callers which found this variant\">'; ",
#      "echo '##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"SAMPLES\">'; ",
#      "echo '##INFO=<ID=SOURCE_POS_END,Number=.,Type=String,Description=\"SOURCE_POS_END\">'; ",
#      "cat /ddn/projects11/got2d/GoT2DSVs/SVGtypes/Dwnldd/t2dgo_chr", chromosome, "_stg1_merged.genotypes.vcf | head -1000 | grep ^#CHROM; ",
#      "grep -v ^# /ddn/projects11/got2d/GoT2DSVs/SVGtypes/Dwnldd/t2dgo_chr", chromosome, "_stg1_merged.genotypes.vcf ;) > ",
#      "/ddn/projects11/got2d/rpearson/SVGtypes/T2DGO/chartl/FinalFreeze/SV.genotyping/t2dgo_chr", chromosome, "_stg1_merged.genotypes.fixed.vcf",
#      sep=""
#    )
#    cat(bashCommand, "\n")
#    system(bashCommand)
    if(!file.exists(sprintf(got2dVcfFmt, chromosome))) {
      cat("Adding extra info lines to vcf\n")
      bashCommand <- paste(
        "(zcat ", sprintf(got2dVcfFmtBeforeFixing, chromosome)," | head -1000 | grep ^##; ",
        "echo '##INFO=<ID=GCFRACTION,Number=1,Type=Float,Description=\"GC content fraction\">'; ",
        "echo '##INFO=<ID=CalledBy,Number=.,Type=String,Description=\"SV callers which found this variant\">'; ",
        "echo '##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"SAMPLES\">'; ",
        "echo '##INFO=<ID=SOURCE_POS_END,Number=.,Type=String,Description=\"SOURCE_POS_END\">'; ",
        "cat ", sprintf(got2dVcfFmtBeforeFixing, chromosome)," | head -1000 | grep ^#CHROM; ",
        "grep -v ^# ", sprintf(got2dVcfFmtBeforeFixing, chromosome)," ;) > ",
        sprintf(got2dVcfFmt, chromosome),
        sep=""
      )
      cat(bashCommand, "\n")
      system(bashCommand)
    }
    
    ############################################################################
    # Load low-pass vcf
    ############################################################################
    cat("Reading low-pass vcf\n")
    vcf <- readVcf(                                                             # read vcf file into a VCF object from package VariantAnnotation
      sprintf(got2dVcfFmt, chromosome),                                         # filename for this chromosome's vcf file
      genome="dummy",                                                           # VCF objects have to include this information, but this is not used in any subsequent analysis, so set to dummy values
      param=ScanVcfParam(geno=c("GT", "FT"))                                    # choose to read in only GT and FT fields from vcf as these are only ones used in subsequent analysis. This keeps RAM usage down
    )
    end(ranges(vcf)) <- values(info(vcf))[["END"]]                              # ranges(vcf) is an object of class GRanges. readVcf populates the start but not the end of the SV, so we have to do it manually. This is important as we use this information later, for example in calculating reciprocal overlaps
    names(rowData(vcf)) <- paste(                                               # here we are making the variant IDs unique by prepending the ID with the chromosome name
      as.character(seqnames(rowData(vcf))),
      names(rowData(vcf)),
      sep="_"
    )
    gc()                                                                        # this is garabge collection to help keep RAM usage down. Probably not necessary but does no harm
    cat("Saving vcf\n")
    save(vcf, file=outputFilenamePreAnnotation)                                 # here we are saving the VCF object in an R object file. This is much quicker to read back into memory than using readVcf to parse the vcf file
#    passVcf <- vcf[values(fixed(vcf))[["FILTER"]] == "PASS"]                    # this would be the way to subset to PASS variant only, but we will do this later on
    
# I had problems when running the following code block so abandoned this
#    ############################################################################
#    # Load thunder vcf
#    ############################################################################
#    cat("Reading thunder vcf\n")
#    thunderVcf <- readVcf(
#      sprintf(thunderVcfFmt, chromosome),
#      genome="dummy",
#      param=ScanVcfParam(fixed="FILTER", info=c("END"))
##      param=ScanVcfParam(fixed="FILTER", info=c("AC", "AN", "END", "GCLENGTH", "GLALTFREQ", "SVLEN", "BAVGPOST", "BRSQ", "LDAF", "AVGPOST", "RSQ"), geno=c("BD"))
#    )
#    end(ranges(thunderVcf)) <- values(info(thunderVcf))[["END"]]
#    gc()
#    cat("Saving vcf\n")
#    save(thunderVcf, file=outputFilenameThunder)
#    
    if(shouldAnnotateWithBatchNum) {
      ############################################################################
      # Load batch1 vcf
      ############################################################################
      cat("Reading batch1 vcf\n")
      batch1vcf <- readVcf(
        sprintf(batch1Fmt, chromosome),
        genome="dummy",
        param=ScanVcfParam(fixed="FILTER", info=c("END"), geno=NA)                # here we choose to read in only the required information from the vcf, e.g. only the END from INFO and none of the FORMAT (genotype) data
      )
      end(ranges(batch1vcf)) <- values(info(batch1vcf))[["END"]]                  
      gc()
      cat("Saving vcf\n")
      save(batch1vcf, file=outputFilenameBatch1)
      
      ############################################################################
      # Load batch2 vcf
      ############################################################################
      cat("Reading batch2 vcf\n")
      batch2vcf <- readVcf(
        sprintf(batch2Fmt, chromosome),
        genome="dummy",
        param=ScanVcfParam(fixed="FILTER", info=c("END"), geno=NA)
      )
      end(ranges(batch2vcf)) <- values(info(batch2vcf))[["END"]]
      gc()
      cat("Saving vcf\n")
      save(batch2vcf, file=outputFilenameBatch2)
      
      ############################################################################
      # Annotate vcf
      ############################################################################
      cat("Annotating vcf\n")
      vcfAnnotated <- annotateSVs(                                                # this function annotates each SV by finding the SV in each of the 4 other vcf files with the highest reciprocal overlap, storing the RO in the INFO, and pulling in all the data from the INFO field for that variant from the other vcf
        vcf,
        list(
  #        thunder                        = thunderVcf,                           # abandoned this due to problems reading Thunder vcfs using readVcf
          batch1                         = batch1vcf,
          batch2                         = batch2vcf,
          thousandGenomesPhase1          = vcf1kg,
          thousandGenomesUsedInDiscovery = vcf1kgUsedInDiscovery
        ),
        annotationStages=c("reciprocalOverlaps")                                  # annotateSVs can do other sorts of annotation, but we're only interested in reciprocal overlaps here
      )
    } else {
      ############################################################################
      # Annotate vcf
      ############################################################################
      cat("Annotating vcf\n")
      vcfAnnotated <- annotateSVs(                                                # this function annotates each SV by finding the SV in each of the 4 other vcf files with the highest reciprocal overlap, storing the RO in the INFO, and pulling in all the data from the INFO field for that variant from the other vcf
        vcf,
        list(
  #        thunder                        = thunderVcf,                           # abandoned this due to problems reading Thunder vcfs using readVcf
#          batch1                         = batch1vcf,
#          batch2                         = batch2vcf,
          thousandGenomesPhase1          = vcf1kg,
          thousandGenomesUsedInDiscovery = vcf1kgUsedInDiscovery
        ),
        annotationStages=c("reciprocalOverlaps")                                  # annotateSVs can do other sorts of annotation, but we're only interested in reciprocal overlaps here
      )
    }

    cat("Fixing annotated vcf\n")
    info(vcfAnnotated) <- values(info(vcfAnnotated))[                           # I have found that VariantAnnotation accumulates columns with names beginning "paramRangeID" which seem to mess up things. I have no idea what these columns are or why they are needed, but these commands seem to remove the associated problems
      ,
      -(grep("paramRangeID", names(values(info(vcfAnnotated)))))
    ]
    fixed(vcfAnnotated) <- values(fixed(vcfAnnotated))[
      ,
      -(grep("paramRangeID", names(values(fixed(vcfAnnotated)))))
    ]
    
#    table(!is.na(GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude)))   # just a debugging command left in for future reference - ignore this
    info(vcfAnnotated) <- BiocGenerics::cbind(                                  # here we annotate SVs by setting a flag in the INFO field for any variant that has any overlap with regions from Adam's file
      values(info(vcfAnnotated)),
      DataFrame(
        "isInAdamsRegions" = !is.na(
          GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude)
        )
      )
    )
#    values(info(vcfAnnotated)) <- cbind(values(info(vcfAnnotated)), DataFrame("isInAdamsRegions" = !is.na(GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude))))   # just a debugging command left in for future reference - ignore this
#    values(info(vcfAnnotated))[["isInAdamsRegions"]] <- rowData(vcfAnnotated) %in% regionsToExclude   # just a debugging command left in for future reference - ignore this
    monomorphicInGoT2d <- which(                                                # here we determine which SVs are monomorphic ref (only non-empty genotype call is "0/0")
      rowSums(
        geno(vcfAnnotated)[["GT"]]=="0/0", na.rm=TRUE
      ) == rowSums(
        !is.na(geno(vcfAnnotated)[["GT"]])
      )
    )
#    table(filt(vcfAnnotated[monomorphicInGoT2d]), useNA="ifany")               # just a debugging command left in for future reference - ignore this
    indexesOfBobHighVPS <- which(                                               # identify which samples are in Bob's list of high variants-per-sample samples
      dimnames(geno(vcfAnnotated)[["GT"]])[[2]] %in% bobHighVPS
    )
    monomorphicInGoT2dAfterRemovingHighVPS <- which(                            # here we determine which SVs are monomorphic ref after high VPS samples have been removed. This is similar to the command a few lines above expect we have removed columns from the matrix corresponding to the high VPS samples
      rowSums(
        geno(vcfAnnotated)[["GT"]][, -(indexesOfBobHighVPS)]=="0/0", na.rm=TRUE
      ) == rowSums(
        !is.na(geno(vcfAnnotated)[["GT"]][, -(indexesOfBobHighVPS)])
      )
    )
    if(length(monomorphicInGoT2dAfterRemovingHighVPS) > 0) {                    # only do the following if there are variants that are monomorphic after removing high VPS samples. This might not be the case for some chromosomes
      filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS]) <- ifelse(     # for variants that are monomorphic after removing high VPS...
        filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS])=="PASS",     # ...change FILTER to NONVARIANT_IN_NONHIGHVPS if currently set to PASS...
        "NONVARIANT_IN_NONHIGHVPS",
        paste(
          filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS]),           # ...or append to current FILTER(s) if not currently set to PASS
          "NONVARIANT_IN_NONHIGHVPS",
          sep=";"
        )
      )
    }
    if(length(which(values(info(vcfAnnotated))[["isInAdamsRegions"]])) > 0) {   # only do the following if there are variants that overlap Adam's regions. This might not be the case for some chromosomes
      filt(                                                                     # change the FILTER only for those variants that we have annotated above as being in Adam's regions
        vcfAnnotated[
          values(info(vcfAnnotated))[["isInAdamsRegions"]]
        ]
      ) <- ifelse(
        filt(
          vcfAnnotated[values(info(vcfAnnotated))[["isInAdamsRegions"]]]
        ) == "PASS",
        "IN_ADAMS_REGIONS",                                                     # ...change FILTER to IN_ADAMS_REGIONS if currently set to PASS...
        paste(                                                                  # ...or append to current FILTER(s) if not currently set to PASS
          filt(vcfAnnotated[values(info(vcfAnnotated))[["isInAdamsRegions"]]]),
          "IN_ADAMS_REGIONS",
          sep=";"
         )
      )
    }
    cat("Saving annotated vcf\n")
    save(vcfAnnotated, file=outputFilenamePostAnnotation)                       # here we are saving the VCF object in an R object file. This is much quicker to read back into memory than using readVcf to parse the vcf file
    
    info(vcfAnnotated) <- values(info(vcfAnnotated))[                           # I have found that VariantAnnotation accumulates columns with names beginning "paramRangeID" which seem to mess up things. I have no idea what these columns are or why they are needed, but these commands seem to remove the associated problems
      ,
      -(grep("paramRangeID", names(values(info(vcfAnnotated)))))
    ]
    fixed(vcfAnnotated) <- values(fixed(vcfAnnotated))[
      ,
      -(grep("paramRangeID", names(values(fixed(vcfAnnotated)))))
    ]
    vcfSites <- VCF(                                                            # here we create a "sites-only" vcf, to keep down file sizes
      rowData = rowData(vcfAnnotated),
      colData = colData(vcfAnnotated),
      exptData = exptData(vcfAnnotated),
      fixed = values(fixed(vcfAnnotated)),
      info = values(info(vcfAnnotated)),
      geno = SimpleList(),
      verbose = FALSE
    )
    writeVcf(vcfSites, filename=outputVcfPostAnnotation, index=TRUE)            # Writes out actual vcf files. The index=TRUE ensures output files and bgzipped and tabix indexed
  } else {
    cat("Loading annotated vcf\n")
    load(outputFilenamePostAnnotation)                                          # if previously created the files, simply load them in rather than doing all the above processing
  }
  gc()                                                                          # this is garabge collection to help keep RAM usage down. Probably not necessary but does no harm
  return(vcfAnnotated)                                                          # return the VCF object to the calling function
}

