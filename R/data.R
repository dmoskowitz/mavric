#' RNA-seq design matrix
#'
#' Design matrix for ENCODE human and mouse RNA-seq data, from GEO accessions GSE36025 and GSE16256, originally published in doi:10.1073/pnas.1413624111, drawn from doi:10.12688/f1000research.6536.1.
#'
#' @docType data
#'
#' @usage data(datasets)
#'
#' @format A data frame with 26 rows (samples) and 4 columns (covariates):
#' \describe{
#'   \item{setname}{Name of sample}
#'   \item{seqBatch}{Sequencing batch in which sample was processed}
#'   \item{species}{Species of sample, either human or mouse}
#'   \item{tissue}{Tissue of origin}
#' }
#'
#' @source \url{http://dx.doi.org/10.5281/zenodo.17606}
"datasets"

#' RNA-seq data matrix
#'
#' Read counts for ENCODE human and mouse RNA-seq data, from GEO accessions GSE36025 and GSE16256, originally published in doi:10.1073/pnas.1413624111, drawn from doi:10.12688/f1000research.6536.1.
#'
#' @docType data
#'
#' @usage data(rawCounts)
#'
#' @format A matrix with 14,744 rows (genes) and 26 columns (samples)
#'
#' @source \url{http://dx.doi.org/10.5281/zenodo.17606}
"rawCounts"
