\name{recoupPlot}
\alias{recoupPlot}
\title{Plot list objects returned by recoup}
\usage{
    recoupPlot(recoupObj, plotParams = list(plot = FALSE, 
        profile = TRUE, heatmap = TRUE, device = "x11"), 
        mainh = 1)
}
\arguments{
    \item{recoupObj}{a list object created from 
    \code{\link{recoup}}.}

    \item{plotParams}{see the \code{plotParams} in the 
    main \code{\link{recoup}} function. A minimum valid
    version is provided for default plotting.}
    
    \item{mainh}{the reference heatmap for ordering
    operations. Normally, calculated in 
    \code{\link{recoup}}. See also the \code{\link{draw}}
    function in the \code{ComplexHeatmap} package. Deafults
    to the first heatmap.}
}
\value{
    This function does not returns anything, just plots the
    \code{recoup} plots.
}
\description{
    This function takes as input argument and output object
    from \code{\link{recoup}} and plots the \code{ggplot2} 
    and \code{ComplexHeatmap} objects stored there.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

pp <- 

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
recoupPlot(test.tss)
}
\author{
    Panagiotis Moulos
}
