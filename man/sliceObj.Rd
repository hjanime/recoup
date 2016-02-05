\name{sliceObj}
\alias{sliceObj}
\title{Subset recoup output list objects}
\usage{
    sliceObj(obj, i = NULL, j = NULL, dropPlots = FALSE,
        rc = NULL) 
}
\arguments{
    \item{obj}{a list object created from 
    \code{\link{recoup}}.}

    \item{i}{vector of numeric or character indices, 
    corresponding to the index or rownames or names of
    reference genomic regions. The \code{design} object
    member will also be subset. If there is a 
    \code{selector} attached to the input object (see
    \code{\link{recoup}} arguments) it will be dropped.}
    
    \item{j}{vector of numeric or character indices
    corresponding to sample index or sample names. These
    will be returned.}
    
    \item{dropPlots}{if profile and/or heatmap plots are
    attached to the input object, they will be recalculated
    if \code{dropPlots=="TRUE"} (default) or dropped 
    otherwise}.
    
    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    A \code{recoup} list object, susbet according to 
    \code{i, j}.
}
\description{
    This function takes as input argument an output object
    from \code{\link{recoup}} and subsets it according to the
    inputs \code{i,j}. The attached plots may or may not be
    recalculated. Other input parameters stores in 
    \code{obj$callopts} are not changed apart from any
    \code{selector} option which is dropped.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

# Calculate coverages
test.tss <- recoup(
    test.input,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    plotParams=list(plot=FALSE,profile=TRUE,
        heatmap=TRUE,device="x11"),
    rc=0.5
)

# Plot coverage profiles
o <- sliceObj(test.tss,i=1:10,j=1)
}
\author{
    Panagiotis Moulos
}
