% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mavric-plotPCs.R
\name{plotPCs}
\alias{plotPCs}
\title{Plot Samples in Subspace of Assigned PCs for a Given Covariate}
\usage{
plotPCs(results, annotation, feature, pcset = NULL, cex = 1,
  lpos = "topleft", transp = TRUE, numberline = T, type = "p")
}
\arguments{
\item{results}{An object output by \code{estVC}}

\item{annotation}{The `annotation' data frame fed to \code{estVC} when creating `results'}

\item{feature}{Column number in `annotation' giving the covariate from which levels and the subspace should be obtained}

\item{pcset}{Rather than using the associated PCs, plot along the specified PCs}

\item{cex}{As in \code{par}, controls the size of points and labels}

\item{lpos}{As in \code{legend}, controls the location of the legend}

\item{transp}{Whether the legend box should be transparent (versus filled)}
}
\value{
None
}
\description{
This function offers a visualization of the separation between levels of a categorical covariate, or across the range of a continuous covariate, within the context of the PC-subspace assigned to the covariate by \code{estVC}.
}
\details{
This function uses the coordinates from the PC-subspace assigned to a particular covariate, and plots the samples, colored according to the value the sample has for that covariate. For categorical covariates, samples are given different color by level; for continuous covariates, samples are plotted in grayscale.

Note, critically, that `feature' refers to the *column number* in `annotation' that is of interest, not to that column's name.
}
\seealso{
\code{estVC} for example usage, \code{plotVars} for Euler diagram plotting, \code{scatterplot3d} in the \code{scatterplot3d} package, \code{par}, and \code{legend}.
}

