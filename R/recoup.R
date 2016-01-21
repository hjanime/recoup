recoup <- function(
    input,
    design=NULL,
    region=c("genebody","tss","tes","custom"),
    type=c("chipseq","rnaseq"),
    genome=c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7","pantro4",
        "susscr3"),
    refdb=c("ensembl","ucsc","refseq"),
    flank=c(2000,2000),
    orderBy=list(
        what=c("none","suma","sumn","maxa","maxn","hcn"),
        order=c("descending","ascending")
    ),
    #orderBy=getDefaultListArgs("orderBy"),
    binParams=list(
        flankBinSize=0,
        regionBinSize=0,
        sumStat=c("mean","median"),
        smooth=TRUE,
        interpolation=c("auto","spline","linear","neighborhood"),
        forceHeatmapBinning=TRUE,
        forcedBinSize=c(50,200)
    ),
    #binParams=getDefaultListArgs("binParams"),
    selector=list(
        id=NULL,
        biotype=getBiotypes(genome),
        exonType=NULL
    ),
    #selector=getDefaultListArgs("selector"),
    preprocessParams=list(
        normalize=c("none","linear","downsample","sampleto"),
        sampleTo=1e+6,
        spliceAction=c("split","keep","remove"),
        spliceRemoveQ=0.75,
        seed=42
    ),
    #preprocessParams=getDefaultListArgs("preprocessParams"),
    plotParams=list(
        plot=TRUE,
        profile=TRUE,
        heatmap=TRUE,
        signalScale=c("natural","log2"),
        heatmapScale=c("common","each"),
        heatmapFactor=1,
        device=c("x11","png","jpg","tiff","bmp","pdf","ps"),
        outputDir=".",
        outputBase=NULL
    ),
    #plotParams=getDefaultListArgs("plotParams"),
    saveParams=list(
        ranges=TRUE,
        coverage=TRUE,
        profile=TRUE,
        profilePlot=TRUE,
        heatmapPlot=TRUE
    ),
    #saveParams=getDefaultListArgs("saveParams"),
    complexHeatmapParams=list(
        main=list(
            cluster_rows=ifelse(length(grep("hc",orderBy$what))>0,TRUE,FALSE),
            cluster_columns=FALSE,
            column_title_gp=grid::gpar(fontsize=12,font=2),
            show_row_names=FALSE,
            show_column_names=FALSE,
            heatmap_legend_param=list(
                color_bar="continuous"
            )
        ),
        group=list(
            cluster_rows=ifelse(length(grep("hc",orderBy$what))>0,TRUE,FALSE),
            cluster_columns=FALSE,
            column_title_gp=grid::gpar(fontsize=12,font=2),
            show_row_names=FALSE,
            show_column_names=FALSE,
            row_title_gp=grid::gpar(fontsize=10,font=2),
            gap=unit(5,"mm"),
            heatmap_legend_param=list(
                color_bar="continuous"
            )
        )
    ),
    #complexHeatmapParams=getDefaultListArgs("complexHeatmapParams"),
    kmParams=list(
        k=0, # Do not perform kmeans
        nstart=20,
        algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen"),
        iterMax=20,
        reference=NULL,
        seed=42
    ),
    #kmParams=getDefaultListArgs("kmParams"),
    strandedParams=list(
        strand=NULL,
        ignoreStrand=TRUE
    ),
    #strandedParams=getDefaultListArgs("strandedParams"),
    bamParams=NULL,
    onTheFly=FALSE, # Directly from BAM w/o Ranges storing, also N/A with BED,
    localDbHome=file.path(path.expand("~"),".recoup"),
    rc=NULL
) {
    if (!is.list(input) && file.exists(input))
        input <- readConfig(input)
    checkInput(input)
    if (is.null(names(input)))
        names(input) <- sapply(input,function(x) return(x$id))
    
    # Check rest of arguments
    region <- tolower(region[1])
    refdb <- tolower(refdb[1])
    type <- tolower(type[1])
    if (!is.null(design) && !is.data.frame(design))
        checkFileArgs("design",design)
    if (!is.data.frame(genome) && file.exists(genome))
        checkFileArgs("genome",genome)
    else if (is.character(genome))
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
    orderByDefault <- getDefaultListArgs("orderBy")
    binParamsDefault <- getDefaultListArgs("binParams")
    selectorDefault <- getDefaultListArgs("selector")
    preprocessParamsDefault <- getDefaultListArgs("preprocessParams")
    plotParamsDefault <- getDefaultListArgs("plotParams")
    saveParamsDefault <- getDefaultListArgs("saveParams")
    kmParamsDefault <- getDefaultListArgs("kmParams")
    strandedParamsDefault <- getDefaultListArgs("strandedParams")
    
    orderBy <- setArg(orderByDefault,orderBy)
    binParams <- setArg(binParamsDefault,binParams)
    if (!is.null(selector))
        selector <- setArg(selectorDefault,selector)
    preprocessParams <- setArg(preprocessParamsDefault,preprocessParams)
    plotParams <- setArg(plotParamsDefault,plotParams)
    saveParams <- setArg(saveParamsDefault,saveParams)
    kmParams <- setArg(kmParamsDefault,kmParams)
    strandedParams <- setArg(strandedParamsDefault,strandedParams)
    
    orderBy <- validateListArgs("orderBy",orderBy)
    binParams <- validateListArgs("binParams",binParams)
    if (!is.null(selector))
        selector <- validateListArgs("selector",selector)
    preprocessParams <- validateListArgs("preprocessParams",preprocessParams)
    plotParams <- validateListArgs("plotParams",plotParams)
    saveParams <- validateListArgs("saveParams",saveParams)
    kmParams <- validateListArgs("kmParams",kmParams)
    strandedParams <- validateListArgs("strandedParams",strandedParams)
    
    if (is.null(plotParams$outputBase))
        plotParams$outputBase <- paste(sapply(input,function(x) return(x$id)),
            collapse="-")
        
    # Check compatibility of orderBy argument and ComplexHeatmap parameters
    # Hierarchical clustering asked in orderBy but otherwise in the heatmap
    # parameters. Clustering is performed.
    if (length(grep("hc",orderBy$what))>0
        && !(complexHeatmapParams$main$cluster_rows 
        || complexHeatmapParams$group$cluster_rows)) {
        warning("Hierarchical clustering asked in the orderBy parameter but ",
            "is set to FALSE in complexHeatmapParams! Will auto-correct to ",
            "perform hierarchical clustering ordering...",immediate.=TRUE)
        complexHeatmapParams$main$cluster_rows <- TRUE
        complexHeatmapParams$group$cluster_rows <- TRUE
    }
    # Hierarchical clustering asked in heatmap parameters but not in the 
    # orderBy directives. Clustering is not performed.
    if ((complexHeatmapParams$main$cluster_rows 
        || complexHeatmapParams$group$cluster_rows)
        && length(grep("hc",orderBy$what))==0) {
        warning("Hierarchical clustering asked in the complexHeatmapParams ",
            "parameter but not in orderBy parameter! Hierarchical clustering ",
            "in the heatmap profile will be turned off",immediate.=TRUE)
        complexHeatmapParams$main$cluster_rows <- FALSE
        complexHeatmapParams$group$cluster_rows <- FALSE
    }
    
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
                    g <- load(file.path(localDbHome,refdb,genome,"gene.rda"))
                    helperRanges <- gene
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
    else {
        rownames(genome) <- as.character(genome[,4])
        genomeRanges <- makeGRangesFromDataFrame(
            df=genome,
            keep.extra.columns=TRUE
        )
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
            design <- tryCatch({
                design[names(genomeRanges),,drop=FALSE]
            },error=function(e) {
                stop("Unexpected error occured! Are you sure that element ",
                    "(row) names in the design file are of the same type as ",
                    "the genome file?")
            },finally={})
        else if (nrow(design)<=length(genomeRanges)) {
            genomeRanges <- tryCatch({
                genomeRanges[rownames(design)]
            },error=function(e) {
                stop("Unexpected error occured! Are you sure that element ",
                    "(row) names in the design file are of the same type as ",
                    "the genome file?")
            },finally={})
            if (type=="rnaseq")
                helperRanges <- tryCatch({
                    helperRanges[rownames(design)]
                },error=function(e) {
                    stop("Unexpected error occured! Are you sure that element ",
                        "(row) names in the design file are of the same type as ",
                        "the genome file?")
                },finally={})
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
        if (type=="chipseq")
            genomeRanges <- applySelectors(genomeRanges,selector,rc=rc)
        if (type=="rnaseq") {
            helperRanges <- applySelectors(helperRanges,selector)
            # Filter and align names if we have helperRanges around
            genomeRanges <- genomeRanges[names(helperRanges)]
        }
    }
    
    # Here we must write code for the reading and normalization of bam files
    # The preprocessRanges function looks if there is a valid (not null) ranges
    # field in input
    #if (!onTheFly)
    input <- preprocessRanges(input,preprocessParams,bamParams=bamParams,rc=rc)

    # Remove unwanted seqnames from reference ranges
    chrs <- unique(unlist(lapply(input,function(x) {
        return(as.character(runValue(seqnames(x$ranges))))
    })))
    if (type=="chipseq") {
        keep <- which(as.character(seqnames(genomeRanges)) %in% chrs)
        genomeRanges <- genomeRanges[keep]
    }
    else if (type=="rnaseq") {
        keeph <- which(as.character(seqnames(helperRanges)) %in% chrs)
        helperRanges <- helperRanges[keeph]
        genomeRanges <- genomeRanges[names(helperRanges)]
        ########################################################################
        ## There must be an R bug with `lengths` here as although it runs in 
        ## Rcmd, it does not pass package building or vignette kniting... But 
        ## for the time being it seems that it is not needed as the name 
        ## filtering works
        #lens <- which(lengths(genomeRanges)==0)
        #if (length(lens)>0)
        #    genomeRanges[lens] <- NULL
        ########################################################################
    }
    
    # Now we must follow two paths according to region type, genebody and custom
    # areas with equal/unequal widths, or tss, tes and 1-width custom areas
    customOne <- FALSE
    if (region=="custom" && all(width(genomeRanges)==1))
        customOne <- FALSE
    if (type=="chipseq")
        input <- coverageRef(input,genomeRanges,region,flank,strandedParams)
            #,bamParams)
    else if (type=="rnaseq")
        input <- coverageRnaRef(input,genomeRanges,helperRanges,flank,
            strandedParams)#,bamParams)
    
    # If normalization method is linear, we must adjust the coverages
    if (preprocessParams$normalize=="linear") {
        linFac <- calcLinearFactors(input,preprocessParams)
        for (n in names(input)) {
            if (linFac[n]==1)
                next
            if (is.null(input[[n]]$coverage$center))
                input[[n]]$coverage <- cmclapply(input[[n]]$coverage,
                    function(x,f) {
                        return(x*f)
                    },linFac[n]
                )
            else {
                input[[n]]$coverage$center <- 
                    cmclapply(input[[n]]$coverage$center,function(x,f) {
                        return(x*f)
                    },linFac[n])
                input[[n]]$coverage$upstream <- 
                    cmclapply(input[[n]]$coverage$upstream,function(x,f) {
                        return(x*f)
                    },linFac[n])
                input[[n]]$coverage$downstream <- 
                    cmclapply(input[[n]]$coverage$downstream,function(x,f) {
                        return(x*f)
                    },linFac[n])
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
            # are normalized (the user is responsible to tell recoup how to do 
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
    
    # Attach some config options for profile and heatmap. irrespectively of 
    # subsequent plotting
    plotOpts <- list(
        xAxisParams=list(
            region=region,
            flank=flank,
            customIsBase=customOne
        ),
        yAxisParams=list(
            signalScale=plotParams$signalScale,
            heatmapScale=plotParams$heatmapScale,
            heatmapFactor=plotParams$heatmapFactor
        ),
        binParams=binParams,
        orderBy=orderBy,
        complexHeatmapParams=complexHeatmapParams
    )
    
    # Coverages and profiles calculated... Now depending on plot option, we go 
    # further or return the enriched input object for saving
    if (!plotParams$profile && !plotParams$heatmap) {
        recoupObj <- toOutput(input,design,saveParams,plotOpts)
        return(recoupObj)
    }
    else
        recoupObj <- toOutput(input,design,
            list(ranges=TRUE,coverage=TRUE,profile=TRUE),plotOpts)
            
    ## Our plot objects
    recoupPlots <- list()
    
    # We must pass the matrices to plotting function
    if (plotParams$profile) {
        message("Constructing genomic coverage profile curve")
        theProfile <- recoupProfile(recoupObj,rc=rc)
        recoupPlots$profilePlot <- theProfile
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
        
        if (binParams$forceHeatmapBinning 
            && (binParams$regionBinSize==0 || binParams$flankBinSize==0)) {
            helpInput <- recoupObj$data
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
            helpInput <- recoupObj$data
        
        helpObj <- recoupObj
        helpObj$data <- helpInput
        message("Constructing genomic coverage heatmap")
        theHeatmap <- recoupHeatmap(helpObj,rc=rc)
        recoupPlots$heatmapPlot <- theHeatmap
        
        # Derive the main heatmap in case of hierarchical clustering        
        mainh <- 1
        if (length(grep("hc",orderBy$what))>0) {
            nc <- nchar(orderBy$what)
            mh <- suppressWarnings(as.numeric(substr(orderBy$what,nc,nc)))
            if (is.na(mh))
                warning("Reference profile for hierarchical clustering order ",
                    "not recognized! Using the 1st...",immediate.=TRUE)
            else if (mh > length(input)) {
                warning("Reference profile (",mh,") for hierarchical ",
                    "clustering order does not exist (the input has only ",
                    length(input)," sources! Using the last...",
                    immediate.=TRUE)
                    mainh <- length(input)
            }
            else
                mainh <- mh
        }
    }
    
    # Overwrite final object so as to return it
    recoupObj <- toOutput(input,design,saveParams,plotOpts,recoupPlots)
    
    # Make make any plots asked
    if (plotParams$plot)
        recoupPlot(recoupObj,plotParams,mainh)
    
    # Return the enriched input object according to save options
    return(recoupObj)
}

coverageRef <- function(input,genomeRanges,region=c("tss","tes","genebody",
    "custom"),flank=c(2000,2000),strandedParams=list(strand=NULL,
    ignoreStrand=TRUE),bamParams=NULL) {
    
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
        input <- coverageAreaRef(input,genomeRanges,region,flank,strandedParams)
            #,bamParams)
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
    bamParams=NULL) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    leftRanges <- getFlankingRanges(mainRanges,flank[1],"upstream")
    rightRanges <- getFlankingRanges(mainRanges,flank[2],"downstream")
    
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        theRanges <- splitBySeqname(input[[n]]$ranges)
        message("Calculating ",region," coverage for ",input[[n]]$name)
        message(" center...")
        input[[n]]$coverage$center <- calcCoverage(theRanges,mainRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" upstream...")
        input[[n]]$coverage$upstream <- calcCoverage(theRanges,leftRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" downstream...")
        input[[n]]$coverage$downstream <- calcCoverage(theRanges,rightRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
    }
    return(input)
}

coverageRnaRef <- function(input,genomeRanges,helperRanges,flank,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),bamParams=NULL) {
    leftRanges <- getFlankingRanges(helperRanges,flank[1],"upstream")
    rightRanges <- getFlankingRanges(helperRanges,flank[2],"downstream")
    
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        theRanges <- splitBySeqname(input[[n]]$ranges)
        message("Calculating genebody coverage for ",input[[n]]$name)
        message(" center...")
        input[[n]]$coverage$center <- calcCoverage(theRanges,genomeRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" upstream...")
        input[[n]]$coverage$upstream <- calcCoverage(theRanges,leftRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand)
        message(" downstream...")
        input[[n]]$coverage$downstream <- calcCoverage(theRanges,rightRanges,
            strand=strandedParams$strand,
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
    if (!is(mask,"GRanges") && !is(mask,"GRangesList"))
        stop("The mask argument must be a GRanges or GRangesList object")
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
    if (is(mask,"GRangesList"))
        x <- mask[[i]]
    else
        x <- mask[i]
    y<-list(
        chromosome=as.character(seqnames(x))[1],
        start=start(x),
        end=end(x),
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
        message(y$chromosome," not found!")
        y$reads <- NULL
    }
    if (length(y$reads)>0) {
        tryCatch({
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            if (length(y$start)>1) { # List of exons, RNA, merge exons
                i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
                    return(s[j]:e[j])
                },y$start,y$end))
                y$coverage <- y$coverage[[cc]][i2k]
            }
            else
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

coverageFromBam <- function(i,mask,input,ignore.strand,pp) {
    if (is(mask,"GRangesList"))
        x <- mask[[i]]
    else
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
    
    switch(pp$spliceAction,
        keep = {
            y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
                param=bp,with.which_label=TRUE),"GRanges")
        },
        split = {
            y$reads <- unlist(grglist(readGAlignments(file=bam.file,
                index=bam.index,param=bp,with.which_label=TRUE)))
        },
        remove = {
            y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
                param=bp,with.which_label=TRUE),"GRanges")
            qu <- quantile(width(y$reads),pp$spliceRemoveQ)
            rem <- which(width(y$reads)>qu)
            if (length(rem)>0)
                y$reads <- y$reads[-rem]
            message("  Excluded ",length(rem)," reads")
        }
    )
    
    if (length(y$reads)>0) {
        tryCatch({
            seqlevels(y$reads) <- as.character(seqnames(x))
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            if (length(y$start)>1) { # List of exons, RNA, merge exons
                i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
                    return(s[j]:e[j])
                },y$start,y$end))
                y$coverage <- y$coverage[[cc]][i2k]
            }
            else
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
                binSize=binParams$regionBinSize,stat=binParams$sumStat,
                interpolation=binParams$interpolation,rc=rc)
            if (binParams$flankBinSize!=0) {
                message(" upstream")
                left <- binCoverageMatrix(input[[n]]$coverage$upstream,
                    binSize=binParams$flankBinSize,stat=binParams$sumStat,
                    interpolation=binParams$interpolation,rc=rc)
                message(" downstream")
                right <- binCoverageMatrix(input[[n]]$coverage$downstream,
                    binSize=binParams$flankBinSize,stat=binParams$sumStat,
                    interpolation=binParams$interpolation,rc=rc)
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
    if (size==0) {
        for (i in 2:length(cvrg)) {
            if (!is.null(cvrg[[i]])) {
                size <- length(cvrg[[i]])
                break
            }
        }
    }
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
    interpolation=c("auto","spline","linear","neighborhood"),rc=NULL) {
    stat <- stat[1]
    interpolation=interpolation[1]
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
    tmp <- cmclapply(tmp,function(x) splitVector(x,binSize,interpolation),rc=rc)
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

applySelectors <- function(ranges,selector,rc=NULL) {
    if (!is.null(selector$id)) {
        ranges <- ranges[selector$ids]
        if (length(ranges)==0)
            stop("No ranges left after using the identifiers provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    if (!is.null(selector$biotype) && !is.null(ranges$biotype)) {
        good <- which(ranges$biotype %in% selector$biotype)
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the biotypes provided with the ",
                "selector field. Are you sure the identifiers between the two ",
                "files are compatible?")
    }
    if (!is.null(selector$exonType) && !is.null(ranges$exon_type)) {
        good <- which(ranges$exonType %in% selector$exonType)
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the exon types provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    return(ranges)
}

toOutput <- function(input,design=NULL,saveParams,plotParams=NULL,
    plotObjs=NULL) {
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
    plots <- list()
    if (!is.null(plotObjs) && saveParams$profilePlot)
        plots$profile <- plotObjs$profilePlot
    if (!is.null(plotObjs) && saveParams$heatmapPlot)
        plots$heatmap <- plotObjs$heatmapPlot
    return(list(
        data=input,
        design=design,
        plotopts=plotParams,
        plots=plots
    ))
}

#applySelectors <- function(ranges,selector,rc=NULL) {
#    if (!is.null(selector$id)) {
#        ranges <- ranges[selector$ids]
#        if (length(ranges)==0)
#            stop("No ranges left after using the identifiers provided with ",
#                "the selector field. Are you sure the identifiers between the ",
#                "two files are compatible?")
#    }
#    if (!is.null(selector$biotype)) {
#        if (is(ranges,"GRanges") && !is.null(ranges$biotype)) {
#            good <- which(ranges$biotype %in% selector$biotype)
#        }
#        else if (is(ranges,"GRangesList") && !is.null(ranges[[1]]$biotype)) {
#            bts <- unlist(lapply(ranges,function(x) {
#                return(x[1]$biotype)
#            }))
#            good <- which(unlist(cmclapply(ranges,function(x) {
#                if (x[1]$biotype %in% type)
#                    return(FALSE)
#                else
#                    return(TRUE)
#            },selector$biotype,rc=rc)))
#        }
#        ranges <- ranges[good]
#        if (length(ranges)==0)
#            stop("No ranges left after using the biotypes provided with the ",
#                "selector field. Are you sure the identifiers between the two ",
#                "files are compatible?")
#    }
#    if (!is.null(selector$exonType) && !is.null(ranges$exon_type)) {
#        if (is(ranges,"GRanges"))
#            good <- which(ranges$exonType %in% selector$exonType)
#        else if (is(ranges,"GRangesList")) {
#            good <- which(unlist(cmclapply(ranges,function(x) {
#                if (x[1]$exon_type %in% type)
#                    return(FALSE)
#                else
#                    return(TRUE)
#            },selector$exonType)))
#        }
#        ranges <- ranges[good]
#        if (length(ranges)==0)
#            stop("No ranges left after using the exon types provided with ",
#                "the selector field. Are you sure the identifiers between the ",
#                "two files are compatible?")
#    }
#    return(ranges)
#}
