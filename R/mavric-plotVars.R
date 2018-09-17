#' Plot Euler Diagram of Variance Structure
#'
#' This function gives a visual estimate of the confounding structure of the data used as input to \code{estVC}. This method uses the assigned PCs to estimate the percentage of variance in the data explained by each covariate, and to what extent there is overlap.
#'
#' This function uses the percent variance explained by each PC (via the eigenvalues) to estimate the percent variance explained by each covariate. Each covariate's variance explained is given by the sum of the percents of its assigned PCs, weighted by the between-cluster sum-of-squares divided by the total sum-of-squares. Overlaps in the Euler diagram represent overlaps in the sets of assigned PCs. Areas are proportional to percentages.
#'
#' Occasionally, the underlying plotting function, \code{venneuler}, has trouble effectively creating the Euler diagram due to complex overlap structures. For this reason, two features have been implemented. First, the method will jitter the underlying values to ensure that no two covariates overlap perfectly. Secondly, the parameter `incvals' controls whether values are added in ascending or descending order. By flipping `incvals,' more attractive outputs are sometimes attainable.
#'
#' @param results An object output by \code{estVC}
#' @param annotation The `annotation' data frame fed to \code{estVC} when creating `results'
#' @param plotzero Whether to include labels for covariates for which no PCs were assigned, and thus for which no variance is explained
#' @param incvals Whether to add covariates to the Euler diagram in ascending order of percent variance explained
#' @export
#' @return None
#' @seealso \code{estVC}, \code{plotPCs} for plotting individual covariates and their assigned subspaces, and \code{venneuler} in the \code{venneuler} package.
plotVars <- function(results, annotation, plotzero = TRUE, incvals = FALSE, plotve = TRUE, maxStress = .1) {
    if(!is.data.frame(annotation)) annotation <- data.frame(annotation)
    results$pcs <- lapply(results$pcs, function(e) e[as.numeric(rownames(e)) > 0, , drop = F])
    vars <- round(unlist(sapply(1:length(results$pcs), function(p) {
        if(nrow(results$pcs[[p]]) == 0) return(0)
        return(sum(results$pcavar[as.numeric(rownames(results$pcs[[p]]))]*sapply(as.numeric(rownames(results$pcs[[p]])), function(i) totVarEx(results$pcasp, annotation, i, p))))
    })))

    vn <- paste(colnames(annotation), paste(vars, "%", sep=''), sep='\n')
    pcex <- sort(unique(unlist(lapply(results$pcs, rownames))))

    pcm <- cbind(matrix(unlist(apply(matrix(unlist(sapply(1:length(results$pcs), function(li) {
        rv <- rep(0, length(pcex))
        if(nrow(results$pcs[[li]]) > 0) {
            e <- rownames(results$pcs[[li]])
            for(i in 1:length(e)) rv[match(e[i], pcex)] <- results$pcavar[as.numeric(e[i])]*totVarEx(results$pcasp, annotation, as.numeric(e[i]), li)
        }
        return(rv)
    })), nrow=length(pcex), dimnames=list(pcex, vn)), 1, function(i) {
            nz <- head(order(i, decreasing=T), n=sum(i != 0))
            if(length(nz) == 1) return(rbind(i[nz], vn[nz]))
            return(rbind(c(diff(t(embed(i[nz], 2))), i[nz[length(nz)]]), sapply(1:length(nz), function(n) paste(vn[nz[1:n]], collapse='&'))))
    })), nrow=2), matrix(sapply(unlist(sapply(1:length(results$pcs), function(li) if(nrow(results$pcs[[li]]) == 0) return(c(0, vn[li])))), function(i) i), nrow=2))

    unacc <- 100-sum(results$pcavar)
    unex <- 100-(unacc+sum(as.numeric(pcm[1,])))
    other.var <- c(unex, unacc)
    names(other.var) <- paste(c("Unexplained", "Discarded"), paste(as.character(round(c(unex, unacc))), "%", sep=''), sep='\n')
    pl.vars <- sort(c(tapply(as.numeric(unlist(pcm[1,])), unlist(lapply(lapply(strsplit(unlist(pcm[2,]), "&"), sort), paste, collapse = '&')), sum), other.var), decreasing = incvals)
    if(!plotzero) pl.vars <- pl.vars[pl.vars > 0]
    if(!plotve) return(pl.vars)
    pl.vars.c <- ceiling(pl.vars)
    repeat {
        pl.tmp <- abs(jitter(pl.vars.c))
        if(plotzero) pl.tmp[pl.vars.c == 0] = 0
        pl.ve <- venneuler::venneuler(pl.tmp)
        if(pl.ve$stress < maxStress) break
    }
    if(plotve) plot(venneuler::venneuler(pl.tmp))
    return(pl.vars)
}
