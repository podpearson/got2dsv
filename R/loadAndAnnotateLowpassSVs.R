# loadAndAnnotateLowpassSVs.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndAnnotateLowpassSVs <- function(
  chromosome                  = 20,
  got2dVcfFmt                 = "/ddn/projects11/got2d/rpearson/SVGtypes/SVG9oct/t2dgo_chr%s_stg1_merged.genotypes.fixed.vcf",
  thunderVcfFmt               = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation/got2d.2874.chr%s.sv.thunder.vcf",
#  thunderVcfFmt               = "/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Dwnldd/INTGfinal/svonly/got2d.2874.chr%s.sv.thunder.vcf.gz",
  batch1Fmt                   = "/ddn/projects11/got2d/GoT2DSVs/discVCFmerge/PASS/1510_chr%s_stg1_6filtered_PASS.discovery.vcf",
  batch2Fmt                   = "/ddn/projects11/got2d/GoT2DSVs/discVCFmerge/PASS/1292_chr%s_stg1_6filtered_PASS.discovery.vcf",
  vcfOmni                     = NULL,
  vcf1kg                      = NULL,
  vcf1kgUsedInDiscovery       = NULL,
  regionsToExclude            = adamsRegionsToExclude("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/Omni_SVraw/Regions_to_drop_082412.txt"),
  bobHighVPS                  = read.delim("/ddn/projects11/got2d/GoT2DSVs/SVGtypes/bobHighLowVariantCountsEmail20121016/high_vps_samples.dat", as.is=TRUE)[["SAMPLE"]],
  outputDir                   = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation",
  outputFilenamePreAnnotation = file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "rda", sep=".")),
  outputFilenamePostAnnotation= file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "annotateSVs.rda", sep=".")),
  outputVcfPostAnnotation     = file.path(outputDir, paste(sprintf(basename(got2dVcfFmt), chromosome), "annotateSVs.vcf", sep=".")),
  outputFilenameThunder       = file.path(outputDir, paste(sprintf(basename(thunderVcfFmt), chromosome), "rda", sep=".")),
  outputFilenameBatch1        = file.path(outputDir, paste(sprintf(basename(batch1Fmt), chromosome), "rda", sep=".")),
  outputFilenameBatch2        = file.path(outputDir, paste(sprintf(basename(batch2Fmt), chromosome), "rda", sep=".")),
  shouldReload                = !file.exists(outputFilenamePostAnnotation)
) {
  if(shouldReload) {
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
    
    ############################################################################
    # Load low-pass vcf
    ############################################################################
    cat("Reading low-pass vcf\n")
    vcf <- readVcf(
      sprintf(got2dVcfFmt, chromosome),
      genome="dummy",
      param=ScanVcfParam(geno=c("GT", "FT"))
    )
    end(ranges(vcf)) <- values(info(vcf))[["END"]]
    gc()
    cat("Saving vcf\n")
    save(vcf, file=outputFilenamePreAnnotation)
#    passVcf <- vcf[values(fixed(vcf))[["FILTER"]] == "PASS"]
    
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
    ############################################################################
    # Load batch1 vcf
    ############################################################################
    cat("Reading batch1 vcf\n")
    batch1vcf <- readVcf(
      sprintf(batch1Fmt, chromosome),
      genome="dummy",
      param=ScanVcfParam(fixed="FILTER", info=c("END"), geno=NA)
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
    vcfAnnotated <- annotateSVs(
      vcf,
      list(
#        thunder                        = thunderVcf,
        batch1                         = batch1vcf,
        batch2                         = batch2vcf,
        thousandGenomesPhase1          = vcf1kg,
        thousandGenomesUsedInDiscovery = vcf1kgUsedInDiscovery
      ),
      annotationStages=c("reciprocalOverlaps")
#      annotationStages=c("setEnds", "reciprocalOverlaps")
    )

    
#    browser()
    cat("Fixing annotated vcf\n")
    info(vcfAnnotated) <- values(info(vcfAnnotated))[, -(grep("paramRangeID", names(values(info(vcfAnnotated)))))]
    fixed(vcfAnnotated) <- values(fixed(vcfAnnotated))[, -(grep("paramRangeID", names(values(fixed(vcfAnnotated)))))]
    
#    table(!is.na(GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude)))
    info(vcfAnnotated) <- BiocGenerics::cbind(values(info(vcfAnnotated)), DataFrame("isInAdamsRegions" = !is.na(GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude)))) # Seems like the bug mentioned above has gone away
#    values(info(vcfAnnotated)) <- cbind(values(info(vcfAnnotated)), DataFrame("isInAdamsRegions" = !is.na(GenomicRanges::match(rowData(vcfAnnotated), regionsToExclude))))
#    values(info(vcfAnnotated))[["isInAdamsRegions"]] <- rowData(vcfAnnotated) %in% regionsToExclude
    monomorphicInGoT2d <- which(rowSums(geno(vcfAnnotated)[["GT"]]=="0/0", na.rm=TRUE) == rowSums(!is.na(geno(vcfAnnotated)[["GT"]])))
    table(filt(vcfAnnotated[monomorphicInGoT2d]), useNA="ifany")
    indexesOfBobHighVPS <- which(dimnames(geno(vcfAnnotated)[["GT"]])[[2]] %in% bobHighVPS)
    monomorphicInGoT2dAfterRemovingHighVPS <- which(rowSums(geno(vcfAnnotated)[["GT"]][, -(indexesOfBobHighVPS)]=="0/0", na.rm=TRUE) == rowSums(!is.na(geno(vcfAnnotated)[["GT"]][, -(indexesOfBobHighVPS)])))
    if(length(monomorphicInGoT2dAfterRemovingHighVPS) > 0) {
      filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS]) <- ifelse(
        filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS])=="PASS",
        "NONVARIANT_IN_NONHIGHVPS",
        paste(filt(vcfAnnotated[monomorphicInGoT2dAfterRemovingHighVPS]), "NONVARIANT_IN_NONHIGHVPS", sep=";")
      )
    }
    if(length(which(values(info(vcfAnnotated))[["isInAdamsRegions"]])) > 0) {
      filt(vcfAnnotated[values(info(vcfAnnotated))[["isInAdamsRegions"]]]) <- ifelse(
        filt(vcfAnnotated[values(info(vcfAnnotated))[["isInAdamsRegions"]]])=="PASS",
        "IN_ADAMS_REGIONS",
        paste(filt(vcfAnnotated[values(info(vcfAnnotated))[["isInAdamsRegions"]]]), "IN_ADAMS_REGIONS", sep=";")
      )
    }
    cat("Saving annotated vcf\n")
    save(vcfAnnotated, file=outputFilenamePostAnnotation)
    
    info(vcfAnnotated) <- values(info(vcfAnnotated))[, -(grep("paramRangeID", names(values(info(vcfAnnotated)))))]
    fixed(vcfAnnotated) <- values(fixed(vcfAnnotated))[, -(grep("paramRangeID", names(values(fixed(vcfAnnotated)))))]
    vcfSites <- VCF(
      rowData = rowData(vcfAnnotated),
      colData = colData(vcfAnnotated),
      exptData = exptData(vcfAnnotated),
      fixed = values(fixed(vcfAnnotated)),
      info = values(info(vcfAnnotated)),
      geno = SimpleList(),
      verbose = FALSE
    )
#    writeVcf(vcfSites, filename=outputVcfPostAnnotation)
    writeVcf(vcfSites, filename=outputVcfPostAnnotation, index=TRUE)
  } else {
#    load(paste("/ddn/projects11/got2d/rpearson/SVGtypes/T2DGO/chartl/FinalFreeze/SV.genotyping/t2dgo_chr", chromosome, "_stg1_merged.genotypes.annotated.vcf.rda", sep=""))
    cat("Loading annotated vcf\n")
    load(outputFilenamePostAnnotation)
  }
  gc()
  return(vcfAnnotated)
}

