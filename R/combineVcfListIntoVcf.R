# combineVcfListIntoVcf.R
# 
# Package: got2dsv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


combineVcfListIntoVcf <- function(
  vcfList,
  infoFieldsToKeep            = names(values(info(vcfList[[1]])))[              # Yes, this is horrible, but it gives the names of the INFO columns that exist in all members of vcfList, and therefore avoids bugs, for example those introduced becuase chr21 vcf has "set" in the INFO field for some unknown reason
    names(values(info(vcfList[[1]]))) %in% names(
      which(
        table(
          unlist(
            lapply(
              vcfList,
              function(vcf) names(values(info(vcf)))
            )
          )
        )==length(
          lapply(
            vcfList,
            function(vcf) names(values(info(vcf)))
          )
        )
      )
    )
  ]
) {
  VCF(
    rowData  = do.call(c, lapply(vcfList, rowData)),
    colData  = colData(vcfList[[1]]),
    exptData = exptData(vcfList[[1]]),
    fixed    = do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) values(fixed(vcf)))),
    info     = do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) values(info(vcf))[, infoFieldsToKeep])),
    geno     = SimpleList(
      sapply(
        names(geno(vcfList[[1]])),
        function(genoField) {
          do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) geno(vcf)[[genoField]]))
        },
        simplify=FALSE,
        USE.NAMES=TRUE
      )
    )
  )
}
