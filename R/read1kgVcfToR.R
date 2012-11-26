# read1kgVcfToR.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


read1kgVcfToR <- function(
  vcfFilename                 = "/ddn/projects11/got2d/GoT2DSVs/1000G_SVgenotypes/ALL.wgs.merged_5_del_call_sets_bps.20101123.sv_dels.low_coverage.sites.vcf",
  outputDir                   = "/ddn/projects11/got2d/rpearson/SVfilteringAndEvaluation",
  outputFilename              = file.path(outputDir, paste(basename(vcfFilename), "rda", sep=".")),
  param                       = ScanVcfParam(fixed=c("FILTER"), geno=NA),
  shouldReload                = !file.exists(outputFilename)
) {
  if(shouldReload) {
    cat("Reading 1kg vcf\n")
    vcf <- readVcf(
      vcfFilename,
      genome="dummy",
      param=param
    )
    end(ranges(vcf)) <- values(info(vcf))[["END"]]
    gc()
    cat("Saving vcf\n")
    save(vcf, file=outputFilename)
  } else {
    cat("Loading 1kg vcf\n")
    load(outputFilename)
  }
  return(vcf)
}
