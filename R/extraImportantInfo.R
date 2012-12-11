# extraImportantInfo.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


extraImportantInfo <- function(
  vcf,
  includesInfoOnBatches       = "batch1_NumberOfOverlaps" %in% names(values(info(vcf))),
  importantInfo               = if(includesInfoOnBatches) {
    c(
      "END", "GCLENGTH", "GLALTFREQ", "GLREFFREQ", "GSELENGTH", "SVLEN", "CalledBy",
      "batch1_NumberOfOverlaps", "batch1_Name", "batch1_Start", "batch1_End", "batch1_Filter", "batch1_ReciprocalOverlap",
      "batch2_NumberOfOverlaps", "batch2_Name", "batch2_Start", "batch2_End", "batch2_Filter", "batch2_ReciprocalOverlap",
      "thousandGenomesPhase1_NumberOfOverlaps", "thousandGenomesPhase1_Name", "thousandGenomesPhase1_Start", "thousandGenomesPhase1_End", "thousandGenomesPhase1_Filter", "thousandGenomesPhase1_ReciprocalOverlap", "EUR_AF",
      "thousandGenomesUsedInDiscovery_NumberOfOverlaps", "thousandGenomesUsedInDiscovery_Name", "thousandGenomesUsedInDiscovery_Start", "thousandGenomesUsedInDiscovery_End", "thousandGenomesUsedInDiscovery_Filter", "thousandGenomesUsedInDiscovery_ReciprocalOverlap", "CalledBy.1",
      "isInAdamsRegions"
    )
  } else {
    c(
      "END", "GCLENGTH", "GLALTFREQ", "GLREFFREQ", "GSELENGTH", "SVLEN", "CalledBy",
      "thousandGenomesPhase1_NumberOfOverlaps", "thousandGenomesPhase1_Name", "thousandGenomesPhase1_Start", "thousandGenomesPhase1_End", "thousandGenomesPhase1_Filter", "thousandGenomesPhase1_ReciprocalOverlap", "EUR_AF",
      "thousandGenomesUsedInDiscovery_NumberOfOverlaps", "thousandGenomesUsedInDiscovery_Name", "thousandGenomesUsedInDiscovery_Start", "thousandGenomesUsedInDiscovery_End", "thousandGenomesUsedInDiscovery_Filter", "thousandGenomesUsedInDiscovery_ReciprocalOverlap", "CalledBy.1",
      "isInAdamsRegions"
    )
  }
) {
  batch <- character(dim(vcf)[1])
  if(includesInfoOnBatches) {
    batch[grep("^b1DEL", names(rowData(vcf)))] <- "GoT2D_batch1"
    batch[grep("^b2DEL", names(rowData(vcf)))] <- "GoT2D_batch2"
  }
  batch[grep("MERGED_DEL", names(rowData(vcf)))] <- "1000genomes"
  
  if(includesInfoOnBatches) {
    newInfoColumns <- DataFrame(
      batch             = batch,
      ROgt0.8withBatch1 = values(info(vcf))[["batch1_ReciprocalOverlap"]]>=0.8,
      ROgt0.8withBatch2 = values(info(vcf))[["batch2_ReciprocalOverlap"]]>=0.8,
      ROgt0.8with1000G  = values(info(vcf))[["thousandGenomesUsedInDiscovery_ReciprocalOverlap"]]>=0.8
    )
  } else {
    newInfoColumns <- DataFrame(
      batch             = batch,
      ROgt0.8with1000G  = values(info(vcf))[["thousandGenomesUsedInDiscovery_ReciprocalOverlap"]]>=0.8
    )
  }
  info(vcf) <- BiocGenerics::cbind(values(info(vcf))[, importantInfo], newInfoColumns)
  vcfSites <- VCF(
    rowData = rowData(vcf),
    colData = colData(vcf),
    exptData = exptData(vcf),
    fixed = values(fixed(vcf)),
    info = values(info(vcf)),
    geno = SimpleList(),
    verbose = FALSE
  )

  if(includesInfoOnBatches) {
    exptData(vcfSites)[["header"]] <- VCFHeader(
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=BiocGenerics::rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="0", Type="Flag", Description="Batch in which the SV was discovered", row.names="batch"),
          DataFrame(Number="0", Type="Flag", Description="Does the SV have a reciprocal overlap >= 0.8 in the GoT2D batch1 discovery set?", row.names="ROgt0.8withBatch1"),
          DataFrame(Number="0", Type="Flag", Description="Does the SV have a reciprocal overlap >= 0.8 in the GoT2D batch2 discovery set?", row.names="ROgt0.8withBatch2"),
          DataFrame(Number="0", Type="Flag", Description="Does the SV have a reciprocal overlap >= 0.8 in the 1000 genomes discovery set?", row.names="ROgt0.8with1000G")
        )
      )
    )
  } else {
    exptData(vcfSites)[["header"]] <- VCFHeader(
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        ALT=fixed(exptData(vcf)[["header"]])[["ALT"]],
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=BiocGenerics::rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="0", Type="Flag", Description="Batch in which the SV was discovered", row.names="batch"),
          DataFrame(Number="0", Type="Flag", Description="Does the SV have a reciprocal overlap >= 0.8 in the 1000 genomes discovery set?", row.names="ROgt0.8with1000G")
        )
      )
    )
  }

  return(vcfSites)
}
