#' Estimate Contributions to Variance
#'
#' This the main end-user function of mavric. Given a data matrix and a set of covariates, this function layers clustering analysis atop PCA to attribute axes of variance to underlying sources.
#'
#' estVC uses a forward selection technique to find a subspace of PCs within which there is separation between levels of categorical covariates.
#'
#' Following PCA, the function iterates over all PCs and calculates both the silhouette statistic, and a null distribution for the silhouette by permuting the PC's values for each sample. Covariates' respective levels are used as cluster assignments when calculating the silhouette. The parameter alpha controls the p-value at which an association between the PC and the covariate is deemed significant. The parameter effsize further controls the attribution by setting a minimum difference between the observed silhouette and the expected value from the permutation distribution. If multiple PCs are more significant than alpha, the PC with the greatest effect size is chosen.
#'
#' The process repeats until no PCs pass the p-value and effect size thresholds. In the first step only, PCs with p-value greater than discthresh are discarded from further steps. In subsequent steps only, only newly added PCs have their values permuted (i.e. not those already selected in previous steps), to test the null hypothesis that adding a given PC provides better separation than would adding a random dimension to the subspace.
#'
#' Given the subspace, the algorithm then forms an axis of variance between pairs, given by the line formed by the pairs' respective centroids in the PC-subspace. Points are then projected onto that axis, and the original measurements are correlated with those projected values to find significant associations.
#'
#' When `autosel' is true, mavric uses a robust LOESS fit to identify measurements with more variance than expected (based on a 95\% confidence interval), given the measurements' respective means.
#'
#' @param data A p x n data matrix, with p measurements on n samples. Typically, this will be count data from an HTS assay.
#' @param annotation An n x l annotation data frame, with l covariates annotated for each of the n samples. This can be a mix of categorical and continuous covariates; categorical covariates are subjected to clustering analysis and continuous covariates are used as the dependent variable in an L1-regularized principal components regression. For categorical covariates, there must be at least two levels, and each level must have at least two samples.
#' @param contrastlist A list describing the contrasts to test. Each element is a list where the first element is the name of covariate, matching a column name in `annotation.' For categorical variables, there should be two additional elements, each with the name of a category (i.e. a level of the factor), indicating that those two categories should be compared; for continuous covariates, the single element of the covariate's name suffices to indicate that testing should occur for high versus low.
#' @param ntop A numeric giving the number of measurements to feed into PCA. For parameter n, the n measurements with the highest variance across all samples will be used.
#' @param sigcor Whether to report those measurements significantly correlated with pairwise differences between categories, for each categorical covariates, and between high/low, for continuous covariates.
#' @param alpha The p-value at which the observed silhouette is deemed sufficient for adding a PC to a subspace. See Details.
#' @param effsize The change in silhouette above the expected value that is deemed sufficient for adding a PC to a subspace. See Details.
#' @param discthresh The p-value threshold above which PCs are discarded from further consideration, during the first forward selection step. See Details.
#' @param ncores The number of CPU cores to use. Note that the memory footprint for each core can be large for large data sets.
#' @param verbose Whether to print out status messages during the run.
#' @param prenorm Whether `data' is already normalized and should not be normalized further.
#' @param autosel Whether to automatically determine which peaks are high-variance, ignoring `ntop.' See Details.
#' @param corMethod Which correlation method to use for identifying significantly correlated measurements. See `method' argument to \code{cor}.
#' @param resampling Whether to use either permutation or bootstrap testing in lieu of the default asymptotic t approximation in the correlation-based statistical testing of the differential analysis.
#' @param pcaobj If mavric has previously been run, the `pcaobj' field of the returned object can be passed here to avoid recomputing it.
#' @param uniqueonly When performing differential comparisons, whether to use all associated PCs, or only those unique to the covariate.
#' @param scaleVar Whether to scale features to unit variance when performing PCA.
#' @param minVar Minimum percent variance a PC must explain to be retained.
#' @param wilcox Whether to use a Wilcoxon rank-sum test rather than MAVRIC's correlation testing.
#' @return A list with the following fields:
#'
#' dat: The `data' object, normalized if applicable.
#'
#' pcaobj: The output of running FactoMineR's PCA on the `ntop' highest-variance measurements of `data' (post-normalization).
#'
#' pcs: A list with a matrix element for each covariate, for the associated PCs and the covariate's mean.
#'
#' pcasp: A matrix giving the coordinates of each sample in PC-space.
#'
#' cor: A list with a list element for each covariate, with each sublist element a matrix giving the correlation between each measurement and the pairwise axis of variance for a given pair of levels in the covariate.
#'
#' @examples
#' ## Human and mouse RNA-seq data from ENCODE, evaluating the amount of
#' ##  variance in gene expression attributable to species-specific
#' ##  versus tissue-specific effects, and quantifying the confounding
#' ##  between species-specific effects and batch effects.
#' ##
#' \dontrun{
#' evc <- estVC(data = rawCounts, annotation = datasets[, -1],
#'  ntop = 10000, sigcor = F, effsize = .01, discthresh = Inf,
#'  ncores = 4)
#' ## Plot variance and confounding amongst experimental covariates
#' plotVars(evc, datasets[, -1])
#' ## Plot subset of PCs best separating samples by species
#' plotPCs(evc, datasets[, -1], "species")
#' }
#'
#' @export
#'
#' @seealso \code{plotPCs} and \code{plotVars} for plotting functionality; \code{silhouette} in the \code{cluster} package for the silhouette statistic; and \code{rlog} and \code{varianceStabilizingTranformation} in the \code{DESeq2} package, and \code{justvsn} in the \code{vsn} package, for variance-stabilizing transformations.
estVC <- function(data, annotation, contrastlist, ntop = Inf, sigcor = TRUE, alpha = .01, effsize = .05, discthresh = 1/3, ncores = 1, verbose = TRUE, prenorm = FALSE, autosel = FALSE, corMethod = c('pearson', 'kendall', 'spearman'), resampling = c('none', 'permutation', 'bootstrap'), pcaobj = NULL, uniqueonly = FALSE, scaleVar = TRUE, minVar = 1, wilcox = FALSE) {
    if(!is.data.frame(annotation)) {
        if(is.factor(annotation) || is.numeric(annotation)) {
            annotation <- data.frame(annotation)
        } else {
            stop("Annotation must be a data frame, factor object, or numeric vector")
        }
    }
    classes <- sapply(1:ncol(annotation), function(j) class(annotation[,j]))
    if(sum(classes == "factor") >= 1) {
        if(verbose) message("Verifying data")
        apply(annotation[, classes == "factor", drop=F], 2, function(j) {
            e <- rle(sort(as.character(j)))$lengths
            if(min(min(e), length(e)) == 1) stop("Each factor must have at least 2 levels, each of which must have at least 2 members")
        })
    }
    corMethod <- match.arg(corMethod)
    resampling <- match.arg(resampling)
    doParallel::registerDoParallel(cores = ncores)
    if(is.null(pcaobj))
        trdv <- htsPCA(data, annotation, ntop, verbose, prenorm, autosel, scaleVar, minVar)
    else
        trdv <- pcaShell(data, pcaobj)
    if(is.vector(trdv$pcasp)) trdv$pcasp <- matrix(trdv$pcasp, ncol=1)
    if(verbose) message("Finding relevant PCs")
    rv <- forSelParRepDisc(trdv$pcasp, annotation, alpha, effsize, discthresh)
    rv$dat <- trdv$dat
    if(!is.null(trdv$as)) rv$as <- trdv$as
    rv$pcavar <- trdv$varex
    rv$pcasp <- trdv$pcasp
    colnames(rv$pcasp) <- sub("^Dim\\.", "PC ", colnames(rv$pcasp))
    rv$varord <- trdv$varord
    rv$pcaobj <- trdv$pcaobj
    names(rv$pcs) <- colnames(annotation)
    if(sigcor) {
        rv$cor <- sigMeasuresProj(rv, annotation, contrastlist, corMethod, resampling, verbose = verbose, uniqueonly = uniqueonly, wilcox = wilcox)
        names(rv$cor) <- colnames(annotation)[colnames(annotation) %in% unlist(lapply(contrastlist, function(e) return(e[1])))]
    }
    return(rv)
}
