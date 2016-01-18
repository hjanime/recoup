\name{preprocessRanges}
\alias{preprocessRanges}
\title{Read and preprocess BAM/BED files to GRanges}
\usage{
    preprocessRanges(input, preprocessParams, 
        bamParams = NULL, rc = NULL)
}
\arguments{
    \item{input}{an input list as in \code{\link{recover}} but
    with the \code{ranges} field of each member filled 
    (e.g. after using \code{\link{preprocessRanges}}).}

    \item{preprocessParams}{see the \code{preprocessParams} 
    argument in the main \code{\link{recover}} function.}
    
    \item{bamParams}{see the \code{bamParams} argument in 
    the main \code{\link{recover}} function.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    This function fills the \code{ranges} field in the
    main \code{input} argument in \code{\link{recover}} 
    function.
}
\description{
    This function reads the BAM/BED files present in the 
    input list object and fills the \code{ranges} field of 
    the latter. At the same time it takes care of certain
    preprocessing steps like normalization.
}
\examples{
\dontrun{
    # This example cannot run as it is dependent on the 
    # presence of BAM files not included with the package
    data("recover_test_data",package="recover")
    pp <- list(normalize="none",sampleTo=1e+6,
        spliceAction="split",seed=42)
    test.input[[1]]$ranges <- NULL
    test.input[[2]]$ranges <- NULL
    test.input <- preprocessRanges(input,pp)
}
}
\author{
    Panagiotis Moulos
}

