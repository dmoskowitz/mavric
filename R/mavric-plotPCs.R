#' Plot Samples in Subspace of Assigned PCs for a Given Covariate
#'
#' This function offers a visualization of the separation between levels of a categorical covariate, or across the range of a continuous covariate, within the context of the PC-subspace assigned to the covariate by \code{estVC}.
#'
#' This function uses the coordinates from the PC-subspace assigned to a particular covariate, and plots the samples, colored according to the value the sample has for that covariate. For categorical covariates, samples are given different color by level; for continuous covariates, samples are plotted in grayscale.
#'
#' Note, critically, that `feature' refers to the *column number* in `annotation' that is of interest, not to that column's name.
#' @param results An object output by \code{estVC}
#' @param annotation The `annotation' data frame fed to \code{estVC} when creating `results'
#' @param feature Column number in `annotation' giving the covariate from which levels and the subspace should be obtained
#' @param pcset Rather than using the associated PCs, plot along the specified PCs
#' @param cex As in \code{par}, controls the size of points and labels
#' @param lpos As in \code{legend}, controls the location of the legend
#' @param transp Whether the legend box should be transparent (versus filled)
#' @export
#' @return None
#' @seealso \code{estVC} for example usage, \code{plotVars} for Euler diagram plotting, \code{scatterplot3d} in the \code{scatterplot3d} package, \code{par}, and \code{legend}.
plotPCs <- function(results, annotation, feature, pcset = NULL, cex = 1, lpos = 'topleft', transp = TRUE, numberline = T, type = 'p') {
    if(transp) transp <- 0
    else transp <- 1
    if(is.null(pcset)) {
        pcset <- as.numeric(rownames(results$pcs[[feature]]))
        pcset <- pcset[pcset >= 1]
        if(length(pcset >= 1) & (class(annotation[,feature]) != "factor")) pcset <- pcset[order(abs(results$pcs[[feature]][as.character(pcset), 1]), decreasing = T)]
    }
    if(length(pcset) == 0) {
        warning("No significantly associated PCs.")
        return(invisible())
    }
    if(length(pcset) > 3) {
        warning("More than three PCs. Only plotting the first three.")
    }
    if(class(annotation[,feature]) == "factor") {
        color.set <- rainbow(length(levels(annotation[, feature])), v = .75)
        color.scheme <- color.set[as.numeric(annotation[, feature])]
    } else {
        color.scheme <- apply(colorRamp(c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"))(1 - (1 * (annotation[, feature] - min(annotation[, feature]))/(max(annotation[, feature]) - min(annotation[, feature])) + 0)), 1, function(i) rgb(i[1], i[2], i[3], maxColorValue = 255))
    }
    if(length(pcset) == 1) {
        par(mar=c(5,5,4,2)+.1)
        if(numberline) {
            mp <- pretty(results$pcasp[, pcset])
            plot(results$pcasp[, pcset], rep(0, nrow(results$pcasp)), axes = F, type = 'n', xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', xlim = range(mp), main = colnames(annotation)[feature], cex.main = cex)
            axis(1, pos = 0, at = mp, cex.axis = cex)
            points(results$pcasp[, pcset], rep(0, nrow(results$pcasp)), col = color.scheme, pch = 16, cex = cex)
            text(max(mp) - diff(range(mp))/2, 0, paste("PC", pcset), pos = 1, offset = 2.5 * par("mgp")[2], cex = cex)
        } else {
            plot(sort(results$pcasp[, pcset]), col=color.scheme[order(results$pcasp[, pcset])], pch=16, main=colnames(annotation)[feature], cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex, ylab = paste("PC", pcset))
        }
    } else if(length(pcset) == 2) {
        axlim <- max(abs(results$pcasp[, pcset]))
        par(mar=c(5,5,4,2)+.1)
        plot(results$pcasp[, pcset], col=color.scheme, pch=16, main=colnames(annotation)[feature], cex = cex, cex.lab = cex, cex.axis = cex, cex.main = cex, cex.main = cex, cex.main = cex, cex.main = cex, xlab = paste("PC", pcset[1]), ylab = paste("PC", pcset[2]))
    } else {
        axlim <- max(abs(results$pcasp[, pcset[1:3]]))
        if(class(annotation[, feature]) == "factor") {
            col <- color.set
            colvar <- as.numeric(annotation[, feature])
            lfeat <- levels(annotation[, feature])
            lfrat <- (length(lfeat) - 1)/length(lfeat)
            colkey <- list(length = .5, cex.axis = .75, at = seq(1 + lfrat/2, length(lfeat) - lfrat/2, by = lfrat), labels = lfeat, dist = -.05)
        } else {
            colvar <- annotation[, feature]
            colkey <- list(length = .5, cex.axis = .75, at = seq(min(annotation[, feature]), max(annotation[, feature]), length.out = 5), dist = -.05)
            col <- rev(colorRampPalette(c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"), space = 'Lab')(100))
        }
        plot3D::scatter3D(results$pcasp[, pcset[1]], results$pcasp[, pcset[2]], results$pcasp[, pcset[3]], colvar = colvar, pch=16, theta = 45, phi = 45, main=colnames(annotation)[feature], cex.lab = cex, cex.axis = cex, cex.symbols = cex, cex.main = cex, xlab = paste("PC", pcset[1]), ylab = paste("PC", pcset[2]), zlab = paste("PC", pcset[3]), nticks = 5, colkey = colkey, ticktype = 'detailed', col = col, type = type)
    }
    if(length(pcset) <= 2) {
        if(class(annotation[,feature]) == "factor") legend(lpos, legend=levels(annotation[, feature]), col = color.set, pch=16, cex = cex, bg=rgb(1,1,1,transp))
        else legend(lpos, legend=c(round(max(annotation[,feature]), 2), "", round(median(annotation[, feature]), 2), "", round(min(annotation[,feature]), 2)), col = apply(colorRamp(c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"))(1 - rev(1 * (quantile(annotation[,feature], (0:4)/4) - min(annotation[,feature]))/(max(annotation[,feature]) - min(annotation[,feature])) + 0)), 1, function(i) rgb(i[1], i[2], i[3], maxColorValue = 255)), lwd=20, cex = cex, bg=rgb(1,1,1,transp))
    }
}
