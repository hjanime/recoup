\name{recoverHeatmap}
\alias{recoverHeatmap}
\title{Plot genomic coverage heatmaps}
\usage{
    recoverHeatmap(recoverObj, rc = NULL)
}
\arguments{
    \item{recoverObj}{a list object created from 
    \code{\link{recover}}.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    The function returns a \code{Heatmap} object with the 
    created heatmap plot. It can be displayed with the 
    \code{draw} function of the \code{ComplexHeatmap}
    package.
}
\description{
    This function takes as input argument and output object
    from \code{\link{recover}} and creates heatmaps 
    depicting genomic coverages using the 
    \code{ComplexHeatmap} package and the options present 
    in the input object. It can be used with saved recover 
    outputs so as to recreate the plots without re-reading 
    BAM/BED files and re-calculating coverages.
}
\examples{
# Load some data
data("recover_test_data",package="recover")

# Calculate coverages
test.tss <- recover(
    test.input,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    plotParams=list(profile=FALSE,heatmap=FALSE),
    rc=0.5
)

# Plot coverage profiles
recoverHeatmap(test.tss,rc=0.5)
}
\author{
    Panagiotis Moulos
}
