splitBySeqname <- function(gr,rc=NULL) {
    #message("Splitting input regions by seqname...")
    gr.list <- cmclapply(levels(seqnames(gr)),function(x,lib) {
        #message("  ",x)
        tmp <- lib[seqnames(lib)==x]
        if (length(tmp)>0) return(tmp) else return(NULL)
    },gr,rc=rc)
    names(gr.list) <- levels(seqnames(gr))
    null <- which(sapply(gr.list,is.null))
    if (length(null)>0)
        gr.list <- gr.list[-null]
    return(gr.list)
}

splitVector <- function(x,n,interp,stat,seed=42) {
    isRle <- ifelse(is(x,"Rle"),TRUE,FALSE)
    if (length(x)<n) {
        if (isRle)
            x <- as.numeric(x)
        switch(interp,
            auto = {
                d <- (n-length(x))/n
                if (d < 0.2) { # Then quite safe for neighborhood method
                    y <- rep(NA,n)
                    set.seed(seed)
                    y[1:2] <- x[1:2]
                    y[(n-1):n] <- x[(length(x)-1):length(x)]
                    orig.pos <- sort(sample(3:(n-2),length(x)-4))
                    y[orig.pos] <- x[3:(length(x)-2)]
                    na <- which(is.na(y))
                    avinds <- lapply(na,function(z) {
                        return(c(z-2,z-1,z+1,z+2))
                    })
                    xx <- unlist(lapply(avinds,function(ii,yy) {
                        return(mean(yy[ii],na.rm=TRUE))
                    },y))
                    y[na] <- xx
                    x <- y
                }
                else { # Spline is the safest
                    x <- spline(x,n=n)$y
                    x[x<0] <- 0
                }
            },
            spline = {
                x <- spline(x,n=n)$y
                x[x<0] <- 0
            },
            inear = {
                x <- approx(x,n=n)$y
                x[x<0] <- 0
            },
            neighborhood = {
                y <- rep(NA,n)
                set.seed(seed)
                y[1:2] <- x[1:2]
                y[(n-1):n] <- x[(length(x)-1):length(x)]
                orig.pos <- sort(sample(3:(n-2),length(x)-4))
                y[orig.pos] <- x[3:(length(x)-2)]
                na <- which(is.na(y))
                avinds <- lapply(na,function(z) {
                    return(c(z-2,z-1,z+1,z+2))
                })
                xx <- unlist(lapply(avinds,function(ii,yy) {
                    return(mean(yy[ii],na.rm=TRUE))
                },y))
                y[na] <- xx
                x <- y
            }
        )
        if (isRle)
            x <- Rle(x)
    }
    bin.size <- floor(length(x)/n)
    dif <- length(x) - bin.size*n 
    bin.fac <- rep(bin.size,n)
    # Random bin increase size to avoid problems
    set.seed(seed)
    add <- sample(1:n,dif)
    bin.fac[add] <- bin.fac[add]+1
    f <- factor(rep(1:n,bin.fac))
    #return(split(x,f))
    S <- split(x,f)
    return(llply(S,stat))
}

readConfig <- function(input) {
    if (missing(input) || !file.exists(input))
        stop("File to read sample info from should be a valid existing text ",
            "file!")
    tab <- read.delim(input)
    if (is.null(tab$id))
        stop("Sample id column not found in ",input,"!")
    if (is.null(tab$file))
        stop("Sample file path column not found in ",input,"!")
    if (is.null(tab$format))
        stop("Sample file format column not found in ",input,"!")
    samples <- as.character(tab$id)
    if (length(samples) != length(unique(samples)))
        stop("Sample identifiers must be unique for each sample!")
    if (length(grep(" ",samples))>0)
        stop("White space is not allowed in sample ids!")
    if (is.null(tab$name))
        nams <- as.character(tab$name)
    else
        nams <- samples
    files <- as.character(tab$file)
    if (any(!file.exists(files))) {
        bi <- which(!file.exists(files))
        stop("Input file ",files[bi]," does not exist! Please check paths...")
    }
    formats <- as.character(tab$format)
    if (!all(formats %in% c("bam","bed")))
        stop("Input formats must be one of \"bam\", \"bed\"")
    cls <- NULL
    if (!is.null(tab$color)) {
        cls <- as.character(tab$color)
        chkcls <- areColors(cls)
        if (!all(chkcls)) {
            warning("Invalid colors found in color column in config file ",
                input,"! Will use automatic colors...")
            cls <- NULL
        }
    }
    output <- vector("list",nrow(tab))
    for (i in 1:nrow(tab)) {
        output[[i]]$id <- samples[i]
        output[[i]]$name <- nams[i]
        output[[i]]$file <- files[i]
        output[[i]]$format <- formats[i]
        output[[i]]$color <- cls[i]
        output[[i]]$ranges <- NULL
        output[[i]]$coverage <- NULL
        output[[i]]$profile <- NULL
    }
    names(output) <- samples
    return(output)
}

kmeansDesign <- function(input,design=NULL,kmParams) {
    if (missing(kmParams)) {
        if (!is.null(input$data))
            kmParams <- getr(input,"kmParams")
        else
            kmParams <- getDefaultListArgs("kmParams")
    }
    else {
        kmParamsDefault <- getDefaultListArgs("kmParams")
        kmParams <- setArg(kmParamsDefault,kmParams)
        kmParams <- validateListArgs("kmParams",kmParams)
    }
    if (!is.null(input$data)) # Fed with recoup object
        input <- input$data
    hasProfile <- sapply(input,function(x) is.null(x$profile))
    if (any(hasProfile))
        stop("Profile matrices for k-means clustering are missing from the ",
            "input object. Have you called the profileMatrix function?")
    if (kmParams$k>0) {
        if (is.null(kmParams$reference)) {
            # Merge matrices to one and perform k-means. As normally coverages 
            # are normalized (the user is responsible to tell recoup how to do 
            # this and has been done at this point), we are legalized to do that
            message("Performing k-means (k=",kmParams$k,") clustering on ",
                "total profiles")
            theBigMatrix <- do.call("cbind",lapply(input,function(x) {
                return(as.matrix(x$profile))
            }))
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
            card <- table(kmorder)[kmorder]
            design$kcluster <- as.factor(paste("Cluster ",kmorder," (",card,
                ")",sep=""))
        }
        else {
            card <- table(kmorder)[kmorder]
            design <- data.frame(kcluster=paste("Cluster ",kmorder," (",card,
                ")",sep=""))
            rownames(design) <- names(kmorder)
        }
    }
    return(design)
}

sliceObj <- function(obj,i=NULL,j=NULL,dropPlots=FALSE,rc=NULL) {
    if (is.null(obj$data))
        stop("No data slot found in obj! Are you sure it's an output from",
            "recoup?")
    if (!is.null(obj$callopts$selector))
        obj$callopts$selector <- NULL
    if (!is.null(i)) {
        if (!is.numeric(i) && !is.character(i))
            stop("Horizontal indexing must be numeric or character!")
        for (s in 1:length(obj$data)) {
            if (!is.null(obj$data[[s]]$coverage))
                obj$data[[s]]$coverage <- obj$data[[s]]$coverage[i]
            if (!is.null(obj$data[[s]]$profile))
                obj$data[[s]]$profile <- obj$data[[s]]$profile[i,,drop=FALSE]
        }
        if (!is.null(obj$design))
            obj$design <- obj$design[i,,drop=FALSE]
    }
    if (!is.null(j)) {
        if (!is.numeric(j) && !is.character(j))
            stop("Vertical indexing must be numeric or character!")
        obj$data <- obj$data[j]
    }
    if (dropPlots)
        obj$plots <- list(profile=NULL,heatmap=NULL)
    else {
        message("Recalculating profile plots after slicing...")
        if (obj$callopts$plotParams$profile) {
            message("  profile")
            obj <- recoupProfile(obj,rc=rc)
        }
        if (obj$callopts$plotParams$heatmap) {
            message("  heatmap")
            obj <- recoupHeatmap(obj,rc=rc)
        }
    }
    return(obj)
}

decideChanges <- function(input,currCall,prevCall) {
    if (is.null(prevCall))
        return(input)
    # Check region and flank
    if (currCall$region != prevCall$region 
        || !all(currCall$flank == prevCall$flank))
        input <- removeData(input,c("coverage","profile"))
    # Check binParams
    if (currCall$binParams$flankBinSize != prevCall$binParams$flankBinSize
        || currCall$binParams$regionBinSize != prevCall$binParams$regionBinSize
        || currCall$binParams$sumStat != prevCall$binParams$sumStat
        || currCall$binParams$smooth != prevCall$binParams$smooth
        || currCall$binParams$interpolation != prevCall$binParams$interpolation
        || currCall$binParams$forceHeatmapBinning != 
            prevCall$binParams$forceHeatmapBinning
        || currCall$binParams$forcedBinSize != prevCall$binParams$forcedBinSize)
        input <- removeData(input,"profile")
    if (currCall$preprocessParams$normalize != 
        prevCall$preprocessParams$normalize
        || currCall$preprocessParams$sampleTo != 
        prevCall$preprocessParams$sampleTo
        || currCall$preprocessParams$spliceAction != 
        prevCall$preprocessParams$spliceAction
        || currCall$preprocessParams$spliceRemoveQ != 
        prevCall$preprocessParams$spliceRemoveQ
        || currCall$preprocessParams$spliceRemoveQ != 
        prevCall$preprocessParams$spliceRemoveQ
        || currCall$preprocessParams$seed != prevCall$preprocessParams$seed)
        input <- removeData(input,c("ranges","coverage","profile"))
        
        return(input)
}

removeData <- function(input,type=c("ranges","coverage","profile")) {
    type <- tolower(type)
    checkTextArgs("type",type,c("ranges","coverage","profile"),multiarg=TRUE)
    if (!is.null(input$data)) # Gave recoup output object
        for (i in 1:length(input$data))
            input$data[[i]][type] <- NULL
    else
        for (i in 1:length(input))
            input[[i]][type] <- NULL
    return(input)
}

calcLinearFactors <- function(input,preprocessParams) {
    hasRanges <- sapply(input,function(x) is.null(x$ranges))
    if (any(hasRanges))
        stop("Please provide input reads before calculation normalization ",
            "factors")
    libSizes <- sapply(input,function(x) return(length(x$ranges)))
    if (preprocessParams$normalize=="linear" 
        || preprocessParams$normalize=="downsample")
        linFac <- min(libSizes)/libSizes
    else if (preprocessParams$normalize=="sampleto")
        linFac <- preprocessParams$sampleTo/libSizes
    names(linFac) <- names(input)
    return(linFac)
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

ssCI <- function(fit) {
    res <- (fit$yin - fit$y)/(1-fit$lev)
    sigma <- sqrt(var(res))
    upper <- fit$y + 3.0*sigma*sqrt(fit$lev)
    lower <- fit$y - 3.0*sigma*sqrt(fit$lev)
    return(list(lower=lower,upper=upper))
}

getDefaultListArgs <- function(arg) {
    switch(arg,
        orderBy = {
            return(list(
                what=c("none","suma","sumn","maxa","maxn","hc"),
                order=c("descending","ascending"),
                custom=NULL
            ))
        },
        binParams = {
            return(list(
                flankBinSize=0,
                regionBinSize=0,
                sumStat=c("mean","median"),
                smooth=TRUE,
                interpolation=c("auto","spline","linear","neighborhood"),
                forceHeatmapBinning=TRUE,
                forcedBinSize=c(50,200)
            ))
        },
        preprocessParams = {
            return(list(
                normalize=c("none","linear","downsample","sampleto"),
                sampleTo=1e+6,
                spliceAction=c("split","keep","remove"),
                spliceRemoveQ=0.75,
                seed=42
            ))
        },
        selector = {
            return(list(
                id=NULL,
                biotype=NULL,
                exonType=NULL
            ))
        },
        strandedParams = {
            return(list(
                strand=NULL,
                ignoreStrand=TRUE
            ))
        },
        plotParams = {
            return(list(
                plot=TRUE,
                profile=TRUE,
                heatmap=TRUE,
                correlation=TRUE,
                signalScale=c("natural","log2"),
                heatmapScale=c("each","common"),
                heatmapFactor=1,
                singleFacet=c("none","wrap","grid"),
                multiFacet=c("wrap","grid"),
                conf=TRUE,
                device=c("x11","png","jpg","tiff","bmp","pdf","ps"),
                outputDir=".",
                outputBase=NULL
            ))
        },
        saveParams = {
            return(list(
                ranges=TRUE,
                coverage=TRUE,
                profile=TRUE,
                profilePlot=TRUE,
                heatmapPlot=TRUE,
                correlationPlot=TRUE
            ))
        },
        kmParams = {
            return(list(
                k=0,
                nstart=20,
                algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen"),
                iterMax=20,
                reference=NULL,
                seed=42
            ))
        }
    )
}

setr <- function(obj,key=c("design","profile","heatmap","correlation","orderBy",
    "kmParams","plotParams"),value) {
    checkTextArgs("key",key,c("design","profile","heatmap","correlation",
        "orderBy","kmParams","plotParams"))
    switch(key,
        design = {
            if (!is.null(design)) {
                if (!is.data.frame(design))
                    design <- read.delim(design,row.names=1)
                nfac <- ncol(design)
                if (length(obj$data)>1 && nfac>2)
                    stop("When more than one files are provided for coverage ",
                        "generation, the maximum number of allowed design ",
                        "factors is 2")
                if (length(obj$data)>1 && nfac>1 && obj$callopts$kmParams$k>0)
                    stop("When more than one files are provided for coverage ",
                        "generation and k-means clustering is also requested, ",
                        "the maximum number of allowed design factors is 1")
                if (length(obj$data)==1 && nfac>3)
                    stop("The maximum number of allowed design factors is 3")
                if (length(obj$data)==1 && nfac>2 && obj$callopts$kmParams$k>0)
                    stop("The maximum number of allowed design factors when ",
                        "k-means clustering is requested is 2")
            }
            obj$design <- design
        },
        profile = {
            if (!all(class(value)==c("gg","ggplot"))) {
                warning("The supplied profile plot is not a ggplot object! ",
                    "Ignoring...",immediate.=TRUE)
                return(obj)
            }
            obj$plots$profile <- value
        },
        heatmap = {
            if (class(value)!=c("HeatmapList")) {
                warning("The supplied heatmap plot is not a Heatmap object! ",
                    "Ignoring...",immediate.=TRUE)
                return(obj)
            }
            obj$plots$heatmap <- value
        },
        correlation = {
            if (!all(class(value)==c("gg","ggplot"))) {
                warning("The supplied correlation plot is not a ggplot ",
                    " object! Ignoring...",immediate.=TRUE)
                return(obj)
            }
            obj$plots$correlation <- value
        },
        orderBy = {
            orderByDefault <- getDefaultListArgs("orderBy")
            orderBy <- setArg(orderByDefault,value)
            orderBy <- validateListArgs("orderBy",orderBy)
            obj$callopts$orderBy <- orderBy
            # Should change also complexHeatmapParams
            if (length(grep("hc",orderBy$what))>0
                && !(obj$callopts$complexHeatmapParams$main$cluster_rows 
                || obj$callopts$complexHeatmapParams$group$cluster_rows)) {
                message("Changing also hierarchical clustering parameters to ",
                    "comply with new orderBy settings.")
                obj$callopts$complexHeatmapParams$main$cluster_rows <- TRUE
                obj$callopts$complexHeatmapParams$group$cluster_rows <- TRUE
            }
            if ((obj$callopts$complexHeatmapParams$main$cluster_rows 
                || obj$callopts$complexHeatmapParams$group$cluster_rows)
                && length(grep("hc",orderBy$what))==0) {
                message("Changing also hierarchical clustering parameters to ",
                    "comply with new orderBy settings.")
                obj$callopts$complexHeatmapParams$main$cluster_rows <- FALSE
                obj$callopts$complexHeatmapParams$group$cluster_rows <- FALSE
            }
        },
        kmParams = {
            kmParamsDefault <- getDefaultListArgs("kmParams")
            kmParams <- setArg(kmParamsDefault,value)
            kmParams <- validateListArgs("kmParams",kmParams)
            obj$callopts$kmParams <- kmParams
        },
        plotParams = {
            plotParamsDefault <- getDefaultListArgs("plotParams")
            plotParams <- setArg(plotParamsDefault,value)
            plotParams <- validateListArgs("plotParams",plotParams)
            obj$callopts$plotParams <- plotParams
        }
    )
    return(obj)
}

getr <- function(obj,key=c("design","profile","heatmap","correlation","orderBy",
    "kmParams","plotParams")) {
    checkTextArgs("key",key,c("design","profile","heatmap","correlation",
        "orderBy","kmParams","plotParams"))
    switch(key,
        design = {
            return(obj$design)
        },
        profile = {
            return(obj$plots$profile)
        },
        heatmap = {
            return(obj$plots$heatmap)
        },
        correlation = {
            return(obj$plots$correlation)
        },
        orderBy = {
            return(obj$callopts$orderBy)
        },
        kmParams = {
            return(obj$callopts$kmParams)
        },
        plotParams = {
            return(obj$callopts$plotParams)
        }
    )
}

getBiotypes <- function(org) {
    if (!(org %in% c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7",
        "pantro4","susscr3")))
        return(NULL)
    switch(org,
        hg18 = {
            return(c("unprocessed_pseudogene","pseudogene","miRNA",
                "retrotransposed","protein_coding","processed_pseudogene",
                "snRNA","snRNA_pseudogene","Mt_tRNA_pseudogene",
                "miRNA_pseudogene","misc_RNA","tRNA_pseudogene","snoRNA",
                "scRNA_pseudogene","rRNA_pseudogene","snoRNA_pseudogene","rRNA","misc_RNA_pseudogene","IG_V_gene","IG_D_gene","IG_J_gene",
                "IG_C_gene","IG_pseudogene","scRNA"))
        },
        hg19 = {
            return(c("pseudogene","lincRNA","protein_coding","antisense",
                "processed_transcript","snRNA","sense_intronic","miRNA",
                "misc_RNA","snoRNA","rRNA","polymorphic_pseudogene",
                "sense_overlapping","3prime_overlapping_ncrna","TR_V_gene",
                "TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene",
                "TR_J_pseudogene","IG_C_gene","IG_C_pseudogene","IG_J_gene",
                "IG_J_pseudogene","IG_D_gene","IG_V_gene","IG_V_pseudogene"))
        },
        hg38 = {
            return(c("protein_coding","polymorphic_pseudogene","lincRNA",
                "unprocessed_pseudogene","processed_pseudogene","antisense",
                "processed_transcript","transcribed_unprocessed_pseudogene",
                "sense_intronic","unitary_pseudogene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","sense_overlapping",
                "transcribed_processed_pseudogene","miRNA","snRNA","misc_RNA",
                "rRNA","snoRNA","IG_J_pseudogene","IG_J_gene","IG_D_gene",
                "3prime_overlapping_ncrna","IG_C_gene","IG_C_pseudogene",
                "pseudogene","TR_V_pseudogene","Mt_tRNA","Mt_rRNA",
                "translated_processed_pseudogene","TR_J_gene","TR_C_gene",
                "TR_D_gene","TR_J_pseudogene","LRG_gene"))
        },
        mm9 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "lincRNA","snoRNA","processed_transcript","misc_RNA","rRNA",
                "sense_overlapping","sense_intronic","polymorphic_pseudogene",
                "non_coding","3prime_overlapping_ncrna","IG_C_gene",
                "IG_J_gene","IG_D_gene","IG_V_gene","ncrna_host"))
        },
        mm10 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "snoRNA","lincRNA","processed_transcript","misc_RNA","rRNA",
                "sense_intronic","sense_overlapping","polymorphic_pseudogene",
                "IG_C_gene","IG_J_gene","IG_D_gene","IG_LV_gene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","TR_V_pseudogene",
                "3prime_overlapping_ncrna"))
        },
        dm3 = {
            return(c("protein_coding","ncRNA","snoRNA","pre_miRNA","pseudogene",
                "snRNA","tRNA","rRNA"))
        },
        rn5 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","misc_RNA"))
        },
        danrer7 = {
            return(c("antisense","protein_coding","miRNA","snoRNA","rRNA",
                "lincRNA","processed_transcript","snRNA","pseudogene",
                "sense_intronic","misc_RNA","polymorphic_pseudogene",
                "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                "non_coding","sense_overlapping"
            ))
        },
        pantro4 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","snRNA","snoRNA","misc_RNA"))
        },
        susscr3 = {
            return(c("antisense","protein_coding","lincRNA","pseudogene",
                "processed_transcript","miRNA","rRNA","snRNA","snoRNA",
                "misc_RNA","non_coding","IG_C_gene","IG_J_gene",
                "IG_V_gene","IG_V_pseudogene"))
        },
        tair10 = {
            return(c("miRNA","ncRNA","protein_coding","pseudogene","rRNA",
                "snoRNA","snRNA","transposable_element","tRNA"))
        }
    )
}

getArg <- function(arg.list,arg.name) {
    return(arg.list[arg.name])
}

setArg <- function(arg.list,arg.name,arg.value=NULL) {
    if (is.list(arg.name))
        arg.list[names(arg.name)] <- arg.name
    else if (is.character(arg.name)) {
        tmp <- vector("list",length(arg.name))
        names(tmp) <- arg.name
        i <- 0
        for (n in arg.name) {
            i <- i + 1
            tmp[[n]] <- arg.value[i]
        }
        arg.list[arg.name] <- tmp
    }
    return(arg.list)
}

areColors <- function(x) {
    return(sapply(x,function(y) {
        tryCatch(is.matrix(col2rgb(y)),
            error=function(e) FALSE)
    }))
}
