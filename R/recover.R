require(parallel)
require(GenomicRanges)
require(GenomicAlignments)
require(rtracklayer)

# input
# region: genebody, tss, tes, custom
# flank: c(-2000,2000)
# design: an optional file with identifiers and split conditions (e.g. high, low)
# genome: hg19, mm9 etc
# 

recover <- function(
    input,
    design=NULL,
    region=c("genebody","tss","tes","custom"),
    type=c("chipseq","rnaseq"),
    genome=c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7","pantro4",
        "susscr3"),
    refdb=c("ensembl","ucsc","refseq"),
    flank=c(2000,2000),
    binParams=list(
        flankBinSize=0,
        regionBinSize=0,
        sumStat=c("mean","median"),
        smooth=TRUE,
        forceHeatmapBinning=TRUE,
        forcedBinSize=c(50,200)
    ),
    selector=list(
        id=NULL,
        biotype=getBiotypes(genome),
        exonType=NULL
    ),
    strandedParams=list(
        strand=NULL,
        ignoreStrand=TRUE
    ),
    preprocessParams=list(
        normalize=c("none","linear","downsample","sampleto"),
        sampleTo=1e+6,
        spliceAction=c("keep","remove","split"),
        seed=42
    ),
    bamParams=NULL,
    plotParams=list(
        profile=TRUE,
        heatmap=TRUE,
        heatmapScale=c("each","common"),
        device=c("x11","png","jpg","tiff","bmp","pdf","ps"),
        outputDir=".",
        outputBase=NULL
    ),
    saveParams=list(
        ranges=TRUE,
        coverage=TRUE,
        profile=FALSE
    ),
    complexHeatmapParams=list(
        main=list(),
        group=list()
    ),
    kmParams=list(
        k=0, # Do not perform kmeans
        nstart=20,
        algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen"),
        iterMax=20,
        reference=NULL,
        seed=42
    ),
    onTheFly=FALSE, # Directly from BAM w/o Ranges storing, also N/A with BED,
    localDbHome=file.path(path.expand("~"),".recover"),
    rc=NULL
) {
    # Check main input
    stored = FALSE
    if (!is.list(input) && file.exists(input)) {
        input <- tryCatch({
            a <- load(input)
            input <- input
            stored = TRUE
        },error=function(e) {
            readConfig(input)
        },finally={})
    }
    checkInput(input)
    names(input) <- sapply(input,function(x) return(x$id))
    
    # Check rest of arguments
    region <- tolower(region[1])
    refdb <- tolower(refdb[1])
    type <- tolower(type[1])
    if (!is.null(design))
        checkFileArgs("design",design)
    if (file.exists(genome))
        checkFileArgs("genome",genome)
    else
        checkTextArgs("genome",genome,c("hg18","hg19","hg38","mm9","mm10","rn5",
            "dm3","danrer7","pantro4","susscr3","tair10"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ensembl","ucsc","refseq"),multiarg=FALSE)
    checkTextArgs("type",type,c("chipseq","rnaseq"),multiarg=FALSE)
    #checkNumArgs("flank",flank,"integer",c(-50000,50000),"both")
    if (any(abs(flank)<500))
        stop("The minimum flanking allowed is 500 bp")
    if (any(abs(flank)>50000))
        stop("The maximum flanking allowed is 50000 bp")
        
    # If type is rnaseq, only genebody plots are valid
    if (type=="rnaseq" && region!="genebody") {
        warning("When type is \"rnaseq\", plots can be created only on ",
            "genebodies! Switching to genebody regions...",immediate.=TRUE)
        region <- "genebody"
    }
    
    # and list arguments
    binParamsDefault <- getDefaultListArgs("binParams")
    preprocessParamsDefault <- getDefaultListArgs("preprocessParams")
    selectorDefault <- getDefaultListArgs("selector")
    strandedParamsDefault <- getDefaultListArgs("strandedParams")
    plotParamsDefault <- getDefaultListArgs("plotParams")
    saveParamsDefault <- getDefaultListArgs("saveParams")
    kmParamsDefault <- getDefaultListArgs("kmParams")
    
    validateListArgs("binParams",binParams)
    validateListArgs("preprocessParams",preprocessParams)
    if (!is.null(selector))
        validateListArgs("selector",selector)
    validateListArgs("strandedParams",strandedParams)
    validateListArgs("plotParams",plotParams)
    validateListArgs("saveParams",saveParams)
    validateListArgs("kmParams",kmParams)
    
    binParams <- setArg(binParamsDefault,binParams)
    preprocessParams <- setArg(preprocessParamsDefault,preprocessParams)
    if (!is.null(selector))
        selector <- setArg(selectorDefault,selector)
    strandedParams <- setArg(strandedParamsDefault,strandedParams)
    plotParams <- setArg(plotParamsDefault,plotParams)
    saveParams <- setArg(saveParamsDefault,saveParams)
    kmParams <- setArg(kmParamsDefault,kmParams)
    
    binParams$sumStat <- tolower(binParams$sumStat[1])
    preprocessParams$normalize <- tolower(preprocessParams$normalize[1])
    preprocessParams$spliceAction <- tolower(preprocessParams$spliceAction[1])
    plotParams$heatmapScale <- tolower(plotParams$heatmapScale[1])
    plotParams$device <- tolower(plotParams$device[1])
    kmParams$algorithm <- kmParams$algorithm[1]
    
    if (is.null(plotParams$outputBase))
        plotParams$outputBase <- paste(sapply(input,function(x) return(x$id)),
            collapse="-")
    checkNumArgs("preprocessParams$sampleTo",preprocessParams$sampleTo,
        "numeric",1e+6,"gte")
    
    if (!is.data.frame(genome)) {
        if (file.exists(genome)) {
            genome <- read.delim(genome) # Must be bed like
            rownames(genome) <- as.character(genome[,4])
            genomeRanges <- makeGRangesFromDataFrame(
                df=genome,
                keep.extra.columns=TRUE
            )
        }
        else {
            # Check if local storage has been set
            if (dir.exists(file.path(localDbHome,refdb,genome))) {
                if (type=="chipseq") {
                    g <- load(file.path(localDbHome,refdb,genome,"gene.rda"))
                    genomeRanges <- gene
                }
                else if (type=="rnaseq") {
                    g <- load(file.path(localDbHome,refdb,genome,
                        "summarized_exon.rda"))
                    genomeRanges <- sexon
                }
            }
            else { # On the fly
                message("Getting annotation on the fly for ",genome," from ",
                    refdb)
                if (type=="chipseq")
                    genome <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
                else if (type=="rnaseq") {
                    helper <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
                    helperRanges <- makeGRangesFromDataFrame(
                        df=helper,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    names(helperRanges) <- as.character(helperRanges$gene_id)
                    genome <- getAnnotation(genome,"exon",refdb=refdb,rc=rc)
                    ann.gr <- makeGRangesFromDataFrame(
                        df=genome,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    message("Merging exons")
                    ann.gr <- reduceExons(ann.gr,multic=multic)
                    names(ann.gr) <- as.character(ann.gr$exon_id)
                    genomeRanges <- split(ann.gr,ann.gr$gene_id)
                }
            }
        }
    }
    
    # Read and check design compatibilities. Check if k-means is requested and
    # message accordingly. If k-means is requested it will be added to the 
    # design data frame
    hasDesign <- FALSE
    if (!is.null(design)) {
        if (!is.data.frame(design))
            design <- read.delim(design,row.names=1)
        nfac <- ncol(design)
        if (length(input)>1 && nfac>2)
            stop("When more than one files are provided for coverage ",
                "generation, the maximum number of allowed design factors is 2")
        if (length(input)>1 && nfac>1 && kmParams$k>0)
            stop("When more than one files are provided for coverage ",
                "generation and k-means clustering is also requested, the ",
                "maximum number of allowed design factors is 1")
        if (length(input)==1 && nfac>3)
            stop("The maximum number of allowed design factors is 3")
        if (length(input)==1 && nfac>2 && kmParams$k>0)
            stop("The maximum number of allowed design factors when k-means ",
                "clustering is requested is 2")
        # Reduce the genomeRanges according to design or the other way
        if (nrow(design)>length(genomeRanges))
            design <- design[names(genomeRanges),,drop=FALSE]
        else if (nrow(design)<=length(genomeRanges)) {
            genomeRanges <- genomeRanges[rownames(design)]
            if (type=="rnaseq")
                helperRanges <- helperRanges[rownames(design)]
        }
        # ...but maybe the names are incompatible
        if (length(genomeRanges)==0)
            stop("No ranges left after using the identifiers provided with ",
                "the design file. Are you sure the identifiers between the ",
                "two files are compatible?")
        if (nrow(design)==0)
            stop("No design elements left after using the identifiers ",
                "provided with the genome file. Are you sure the identifiers ",
                "between the two files are compatible?")
    }
    
    # Apply the rest of the filters if any to reduce later computational burden
    if (!is.null(selector)) {
        genomeRanges <- applySelectors(genomeRanges,selector)
        if (type=="rnaseq")
            helperRanges <- applySelectors(helperRanges,selector)
    }
    
    # Align names if we have helperRanges around
    helperRanges <- helperRanges[names(genomeRanges)]
    
    # Here we must write code for the reading and normalization of bam files
    # The preprocessRanges function looks if there is a valid (not null) ranges
    # field in input
    if (!onTheFly)
        input <- preprocessRanges(input,preprocessParams,bamParams=bamParams,
            rc=rc)
    
    # Now we must follow two paths according to region type, genebody and custom
    # areas with equal/unequal widths, or tss, tes and 1-width custom areas
    customOne <- FALSE
    if (region=="custom" && all(width(genomeRanges)==1))
        customOne <- FALSE
    input <- coverageRef(input,genomeRanges,region,flank,strandedParams)#,bamParams)
    
    # If normalization method is linear, we must adjust the coverages
    if (preprocessParams$normalize == "linear") {
        linFac <- calcLinearFactors(input,preprocessParams)
        for (n in names(input)) {
            if (is.null(input[[n]]$coverage$center))
                input[[n]]$coverage <- input[[n]]$coverage * linFac[n]
            else {
                input[[n]]$coverage$center <- 
                    input[[n]]$coverage$center*linFac[n]
                input[[n]]$coverage$upstream <- 
                    input[[n]]$coverage$upstream*linFac[n]
                input[[n]]$coverage$downstream <- 
                    input[[n]]$coverage$downstream*linFac[n]
            }
        }
    }
    
    # Now we must summarize and create the matrices. If genebody or unequal 
    # custom lengths, bin is a must, else we leave to user
    mustBin <- FALSE
    if (region=="genebody")
        mustBin <- TRUE
    if (region=="custom") {
        w <- width(genomeRanges)
        if (any(w!=w[1]))
            mustBin <- TRUE
    }
    
    if (mustBin) {
        if (binParams$regionBinSize==0) {
            warning("Central region bin size not set for a region that must ",
                "be binned! Setting to 1000...",immediate.=TRUE)
            binParams$regionBinSize <- 1000
        }
    }
    input <- profileMatrix(input,region,binParams,rc)
    
    # Perform the k-means clustering if requested and append to design (which
    # has been checked, if we are allowed to do so)
    if (kmParams$k>0) {
        if (is.null(kmParams$reference)) {
            # Merge matrices to one and perform k-means. As normally coverages 
            # are normalized (the user is responsible to tell recover how to do 
            # this and has been done at this point), we are legalized to do that
            message("Performing k-means (k=",kmParams$k,") clustering on ",
                "total profiles")
            theBigMatrix <- do.call("cbind",lapply(input,function(x) {
                return(x$profile)
            }))
            set.seed(kmParams$seed)
            kcl <- kmeans(theBigMatrix,centers=kmParams$k,
                iter.max=kmParams$iterMax,nstart=kmParams$nstart,
                algorithm=kmParams$algorithm)
        }
        else {
            message("Performing k-means (k=",kmParams$k,") clustering using ",
                "the ",input[[kmParams$reference]]$name," sample profile as ",
                "reference")
            set.seed(kmParams$seed)
            kcl <- kmeans(input[[kmParams$reference]]$profile,
                centers=kmParams$k,iter.max=kmParams$iterMax,
                nstart=kmParams$nstart,algorithm=kmParams$algorithm)
        }
        kmorder <- kcl$cluster
        #names(kmorder) <- rownames(input[[1]]$profile)
        if (!is.null(design)) {
            kmorder <- kmorder[rownames(design)]
            design$kcluster <- as.factor(paste("Cluster",kmorder))
        }
        else {
            design <- data.frame(kcluster=paste("Cluster",kmorder))
            rownames(design) <- names(kmorder)
        }
    }
    
    # Coverages and profiles calculated... Now depending on plot option, we go 
    # further or return the enriched input object for saving
    if (!plotParams$profile && !plotParams$heatmap)
        return(toOutput(input,saveParams))
    
    # Attach some config options for profile and heatmap
    plotOpts <- list(
        region=region,
        type=type,
        flank=flank,
        binParams=binParams,
        customOne=customOne,
        heatmapScale=plotParams$heatmapScale,
        complexHeatmapParams=complexHeatmapParams
    )
    for (n in names(input))
        input[[n]]$plotopts <- plotOpts
    
    # We must pass the matrices to plotting function
    if (plotParams$profile) {
        message("Constructing genomic coverage profile curve")
        theProfile <- recoverProfile(input,design=design,rc=rc)
        if (plotParams$device=="x11")
            plot(theProfile)
        else
            ggsave(filename=paste(plotParams$outputBase,"_profile.",
                plotParams$device,sep=""),plot=theProfile,
                path=plotParams$outputDir)
    }
    
    # Some default binning MUST be applied for the heatmap... Otherwise it could
    # take insanely long time and space to draw/store
    if (plotParams$heatmap) {
        # Inform the user about enforced binning (or not)
        if (region %in% c("tss","tes") || customOne) {
            if (binParams$regionBinSize==0 && binParams$forceHeatmapBinning) {
                message("The resolution of the requested profiles will be ",
                    "lowered to avoid increased computation time and/or ",
                    "storage space for heatmap profiles...")
                
            }
            else if (binParams$regionBinSize==0
                && !binParams$forceHeatmapBinning)
                warning("forceHeatmapBinning is turned off for high ",
                    "resolution plotting. Be prepared for long computational ",
                    "times and big figures!",immediate.=TRUE)
        }
        else {
            if ((binParams$regionBinSize==0 || binParams$flankBinSize==0)
                && binParams$forceHeatmapBinning) {
                message("The resolution of the requested profiles will be ",
                    "lowered to avoid increased computation time and/or ",
                    "storage space for heatmap profiles...")
                
            }
            else if ((binParams$regionBinSize==0 || binParams$flankBinSize==0)
                && !binParams$forceHeatmapBinning)
                warning("forceHeatmapBinning is turned off for high ",
                    "resolution plotting. Be prepared for long computational ",
                    "times and big figures!",immediate.=TRUE)
        }
        
        if (binParams$forceHeatmapBinning) {
            helpInput <- input
            if (region %in% c("tss","tes") || customOne) {
                for (n in names(helpInput)) {
                    message("Calculating ",region," profile for ",
                        helpInput[[n]]$name)
                    helpInput[[n]]$profile <- 
                        binCoverageMatrix(helpInput[[n]]$coverage,
                            binSize=binParams$forcedBinSize[2],
                            stat=binParams$sumStat,rc=rc)
                }
            }
            else {
                for (n in names(helpInput)) {
                    message("Calculating ",region," profile for ",
                        helpInput[[n]]$name)
                    message(" center")
                    center <- binCoverageMatrix(helpInput[[n]]$coverage$center,
                        binSize=binParams$forcedBinSize[2],
                        stat=binParams$sumStat,rc=rc)                   
                    message(" upstream")
                    left <- binCoverageMatrix(helpInput[[n]]$coverage$upstream,
                        binSize=binParams$forcedBinSize[1],
                        stat=binParams$sumStat,rc=rc)
                    message(" downstream")
                    right <- binCoverageMatrix(
                        helpInput[[n]]$coverage$downstream,
                        binSize=binParams$forcedBinSize[1],
                        stat=binParams$sumStat,rc=rc)
                    helpInput[[n]]$profile <- cbind(left,center,right)
                }
            }
        }
        else
            helpInput <- input
        
        message("Constructing genomic coverage heatmap")
        theHeatmap <- recoverHeatmap(helpInput,design=design,rc=rc)
        
        if (plotParams$device=="x11") {
            dev.new()
            draw(theHeatmap,gap=unit(1,"cm"))
        }
        else {
            # Starting from width=4, we add 1.5 inches for each heatmap
            iw <- 4 + (length(input)-1)*1.5
            if (plotParams$device == "pdf")
                graphicsOpen(plotParams$device,paste(plotParams$outputBase,
                    "_heatmap.",plotParams$device,sep=""),width=iw)
            else
                graphicsOpen(plotParams$device,paste(plotParams$outputBase,
                    "_heatmap.",plotParams$device,sep=""),width=iw,units="in")
            draw(theHeatmap,gap=unit(1,"cm"))
            graphicsClose(plotParams$device)
        }
    }
    
    # Return the enriched input object according to save options
    return(toOutput(input,saveParams))
}

coverageRef <- function(input,genomeRanges,region=c("tss","tes","genebody",
    "custom"),flank=c(2000,2000),strandedParams=list(strand=NULL,
    ignoreStrand=TRUE),bamParams=TRUE) {
    
    if (region %in% c("tss","tes","custom")) {
        if (region %in% c("tss","tes"))
            input <- coverageBaseRef(input,genomeRanges,region,flank,
                strandedParams)#,bamParams)
        else if (region=="custom") {
            if (all(width(genomeRanges)==1))
                input <- coverageBaseRef(input,genomeRanges,region,flank,
                    strandedParams)#,bamParams)
            else
                input <- coverageAreaRef(input,genomeRanges,region,flank,
                    strandedParams)#,bamParams)
        }
    }
    else # is genebody
        input <- coverageAreaRef(input,genomeRanges,region,flank,strandedParams)#,
            #bamParams)
    return(input)
}

coverageBaseRef <- function(input,genomeRanges,region,flank,strandedParams) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating ",region," coverage for ",input[[n]]$name)
        if (!is.null(input[[n]]$ranges))
            input[[n]]$coverage <- calcCoverage(input[[n]]$ranges,mainRanges,
                strand=strandedParams$strand,
                ignore.strand=strandedParams$ignoreStrand)
        else
            input[[n]]$coverage <- calcCoverage(input[[n]]$file,mainRanges,
                strand=strandedParams$strand,
                ignore.strand=strandedParams$ignoreStrand)
    }
    return(input)
}

coverageAreaRef <- function(input,genomeRanges,region,flank,strandedParams,
    bamParams) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    leftRanges <- getFlankingRanges(mainRanges,flank[1],"upstream")
    rightRanges <- getFlankingRanges(mainRanges,flank[2],"downstream")
    
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating ",region," coverage for ",input[[n]]$name)
        message(" center...")
        input[[n]]$coverage$center <- calcCoverage(input[[n]]$ranges,mainRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" upstream...")
        input[[n]]$coverage$upstream <- calcCoverage(input[[n]]$ranges,
            leftRanges,strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" downstream...")
        input[[n]]$coverage$downstream <- calcCoverage(input[[n]]$ranges,
            rightRanges,strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
    }
    return(input)
}

coverageRnaRef <- function(input,genomeRanges,helperRanges,flank,strandedParams,
    bamParams) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    leftRanges <- getFlankingRanges(helperRanges,flank[1],"upstream")
    rightRanges <- getFlankingRanges(helperRanges,flank[2],"downstream")
    
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating ",region," coverage for ",input[[n]]$name)
        message(" center...")
        input[[n]]$coverage$center <- calcCoverage(input[[n]]$ranges,mainRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" upstream...")
        input[[n]]$coverage$upstream <- calcCoverage(input[[n]]$ranges,
            leftRanges,strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" downstream...")
        input[[n]]$coverage$downstream <- calcCoverage(input[[n]]$ranges,
            rightRanges,strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
    }
    return(input)
}

getRegionalRanges <- function(ranges,region,flank) {
    switch(region,
        genebody = {
            return(ranges)
        },
        tss = {
            return(promoters(ranges,upstream=flank[1],downstream=flank[2]))
        },
        tes = {
            tmp <- resize(ranges,width=1,fix="end")
            return(promoters(tmp,upstream=flank[1],downstream=flank[2]))
        },
        custom = {
            if (all(width(ranges)==1))
                return(promoters(ranges,upstream=flank[1],downstream=flank[2]))
            else
                return(ranges)
        }
    )
}

getFlankingRanges <- function(ranges,flank,dir=c("upstream","downstream")) {
    dir = dir[1]
    if (dir=="upstream")
        return(promoters(ranges,upstream=flank,downstream=0))
    else if (dir=="downstream") {
        return(flank(ranges,width=flank,start=FALSE,both=FALSE))
    }
}

calcCoverage <- function(input,mask,strand=NULL,ignore.strand=TRUE) {
    if (!is(input,"GRanges") && !is.list(input) && is.character(input)
        && !file.exists(input))
        stop("The input argument must be a GenomicRanges object or a valid ",
            "BAM file or a list of GenomicRanges")
    if (!is(mask,"GRanges"))
        stop("The mask argument must be a GenomicRanges object")
    isBam <- FALSE
    if (is.character(input) && file.exists(input))
        isBam <- TRUE
    if (!is.null(strand) && !is.list(strand) && !isBam) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    if (!is.list(input) && !isBam)
        input <- splitBySeqname(input)
    index <- 1:length(mask)
    #message("Calculating coverage...")
    if (isBam)
        cov <- cmclapply(index,coverageFromBam,mask,input,ignore.strand)
    else
        cov <- cmclapply(index,coverageFromRanges,mask,input,ignore.strand)
    names(cov) <- names(mask)
    gc(verbose=FALSE)
    #message("Done!")
    return(cov) # Rle
}

coverageFromRanges <- function(i,mask,input,ignore.strand) {
    x <- mask[i]
    y<-list(
        chromosome=as.character(seqnames(x)),
        start=start(x),
        end=end(x),
        strand=as.character(strand(x)),
        reads=NULL,
        coverage=NULL
    )
    if (!is.null(input[[y$chromosome]])) {
        y$reads <- input[[y$chromosome]][
            subjectHits(findOverlaps(x,input[[y$chromosome]],
                ignore.strand=ignore.strand))]
    }
    else {
        message(y$chromosome,"not found!")
        y$reads <- NULL
    }
    if (length(y$reads)>0) {
        tryCatch({
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        },
        error=function(e) {
            message("Caught invalid genomic area!")
            print(mask[i])
            message("Will return zero coverage")
            return(NULL)
        },finally={})
    }
    else
        return(NULL)
}

coverageFromBam <- function(i,mask,input,ignore.strand) {
    x <- mask[i]
    y<-list(
        chromosome=as.character(seqnames(x)),
        start=start(x),
        end=end(x),
        strand=as.character(strand(x)),
        reads=NULL,
        coverage=NULL
    )
    bam.file <- input
    bam.index <- paste(input,"bai",sep=".")
    bp <- ScanBamParam(which=x)
    y$reads <- as(readGAlignments(file=bam.file,index=bam.index,param=bp,
        with.which_label=TRUE),"GRanges")
    
    if (length(y$reads)>0) {
        tryCatch({
            seqlevels(y$reads) <- as.character(seqnames(x))
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        },
        error=function(e) {
            message("Caught invalid genomic area!")
            print(mask[i])
            message("Will return zero coverage")
            return(NULL)
        },finally={})
    }
    else
        return(NULL)
}

coverageFromRnaRanges <- function(i,mask,input,ignore.strand) {
    x <- mask[i] # Is a set of ranges, not just one
    y<-list(
        chromosome=as.character(seqnames(x))[1],
        #start=start(x),
        #end=end(x),
        strand=as.character(strand(x))[1],
        reads=NULL,
        coverage=NULL
    )
    if (!is.null(input[[y$chromosome]])) {
        y$reads <- input[[y$chromosome]][
            subjectHits(findOverlaps(x,input[[y$chromosome]],
                ignore.strand=ignore.strand))]
    }
    else {
        message(y$chromosome,"not found!")
        y$reads <- NULL
    }
    if (length(y$reads)>0) {
        tryCatch({
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        },
        error=function(e) {
            message("Caught invalid genomic area!")
            print(mask[i])
            message("Will return zero coverage")
            return(NULL)
        },finally={})
    }
    else
        return(NULL)
}

profileMatrix <- function(input,region,binParams,rc=NULL) {
    if (!is.null(input[[1]]$coverage$center)) {
        for (n in names(input)) {
            message("Calculating ",region," profile for ",input[[n]]$name)
            message(" center")
            center <- binCoverageMatrix(input[[n]]$coverage$center,
                binSize=binParams$regionBinSize,stat=binParams$sumStat,rc=rc)
            if (binParams$flankBinSize!=0) {
                message(" upstream")
                left <- binCoverageMatrix(input[[n]]$coverage$upstream,
                    binSize=binParams$flankBinSize,stat=binParams$sumStat,rc=rc)
                message(" downstream")
                right <- binCoverageMatrix(input[[n]]$coverage$downstream,
                    binSize=binParams$flankBinSize,stat=binParams$sumStat,rc=rc)
            }
            else {
                message(" upstream")
                left <- baseCoverageMatrix(input[[n]]$coverage$upstream,rc=rc)
                message(" downstream")
                right <- baseCoverageMatrix(input[[n]]$coverage$downstream,
                    rc=rc)
            }
            input[[n]]$profile <- cbind(left,center,right)
        }
    }
    else {
        for (n in names(input)) {
            message("Calculating ",region," profile for ",input[[n]]$name)
            if (binParams$regionBinSize!=0)
                input[[n]]$profile <- 
                    binCoverageMatrix(input[[n]]$coverage,
                        binSize=binParams$regionBinSize,stat=binParams$sumStat,
                            rc=rc)
            else
                input[[n]]$profile <- 
                    baseCoverageMatrix(input[[n]]$coverage,rc=rc)
        }
    }
    return(input)
}

baseCoverageMatrix <- function(cvrg,rc=NULL) {
    size <- length(cvrg[[1]])
    tmp <- cmclapply(cvrg,function(x) {
        if (class(x)=="Rle")
            return(as.numeric(x))
    },rc=rc)
    null <- which(sapply(tmp,is.null))
    if (length(null)>0) {
        for (j in null) {
            fill <- rep(0,size)
            tmp[[j]] <- fill
        }
    }
    return(do.call("rbind",tmp))
}

binCoverageMatrix <- function(cvrg,binSize=1000,stat=c("mean","median"),
    rc=NULL) {
    stat <- stat[1]
    tmp <- cmclapply(cvrg,function(x) {
        if (class(x)=="Rle")
            return(as.numeric(x))
    },rc=rc)
    null <- which(sapply(tmp,is.null))
    if (length(null)>0) {
        for (j in null) {
            fill <- rep(0,binSize)
            tmp[[j]] <- fill
        }
    }
    tmp <- cmclapply(tmp,function(x) splitVector(x,binSize),rc=rc)
    if (stat=="mean") {
        statMatrix <- do.call("rbind",cmclapply(tmp,function(x) 
            unlist(lapply(x,function(y) mean(y))),rc=rc))
    }
    else if (stat=="median") {
        statMatrix <- do.call("rbind",cmclapply(tmp,function(x) 
            unlist(lapply(x,function(y) median(y))),rc=rc))
    }
    return(statMatrix)
}

applySelectors <- function(ranges,selector) {
    if (!is.null(selector$id)) {
        ranges <- ranges[selector$ids]
        if (length(ranges)==0)
            stop("No ranges left after using the identifiers provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    if (!is.null(selector$biotype) && !is.null(ranges$biotype)) {
        if (is(ranges,"GRanges"))
            good <- which(ranges$biotype %in% selector$biotype)
        else if (is(ranges,"GRangesList")) {
            good <- which(unlist(cmclapply(ranges,function(x) {
                if (x[1]$biotype %in% type)
                    return(FALSE)
                else
                    return(TRUE)
            },selector$biotype)))
        }
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the biotypes provided with the ",
                "selector field. Are you sure the identifiers between the two ",
                "files are compatible?")
    }
    if (!is.null(selector$exonType) && !is.null(ranges$exon_type)) {
        if (is(ranges,"GRanges"))
            good <- which(ranges$exonType %in% selector$exonType)
        else if (is(ranges,"GRangesList")) {
            good <- which(unlist(cmclapply(ranges,function(x) {
                if (x[1]$exon_type %in% type)
                    return(FALSE)
                else
                    return(TRUE)
            },selector$exonType)))
        }
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the exon types provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    return(ranges)
}

toOutput <- function(input,saveParams) {
    if (!saveParams$ranges) {
        for (i in 1:length(input))
            input[[i]]$ranges <- NULL
    }
    if (!saveParams$coverage) {
        for (i in 1:length(input))
            input[[i]]$coverage <- NULL
    }
    if (!saveParams$profile) {
        for (i in 1:length(input))
            input[[i]]$profile <- NULL
    }
    return(input)
}
