test_recover <- function() {
    data("recover_test_data",package="recover")

    test.tss <- recover(
        test.input,
        design=NULL,
        region="tss",
        type="chipseq",
        genome=test.genome,
        flank=c(2000,2000),
        selector=NULL,
        rc=0.5
    )
    
    test.gb <- recover(
        test.input,
        design=test.design,
        region="genebody",
        type="chipseq",
        genome=test.genome,
        flank=c(2000,2000),
        binParams=list(flankBinSize=50,regionBinSize=150),
        orderBy=list(what="hc1"),
        selector=NULL,
        rc=0.5
    )
    
    checkTrue(is.list(test.tss))
    checkTrue(!is.null(test.tss$data))
    checkTrue(is.list(test.gb))
    checkTrue(!is.null(test.gb$data))
}
