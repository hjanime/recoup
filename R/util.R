splitBySeqname <- function(gr) {
    #message("Splitting input regions by seqname...")
    gr.list <- cmclapply(levels(seqnames(gr)),function(x,lib) {
        #message("  ",x)
        tmp <- lib[seqnames(lib)==x]
        if (length(tmp)>0) return(tmp) else return(NULL)
    },gr)
    names(gr.list) <- levels(seqnames(gr))
    null <- which(sapply(gr.list,is.null))
    if (length(null)>0)
        gr.list <- gr.list[-null]
    return(gr.list)
}

splitVector <- function(x,n,seed=42) {
    # If the length of the original vector is smaller, randomly distribute some
    # zeros along x and issue a warning
    if (length(x)<n) {
        y <- rep(0,n)
        set.seed(seed)
        orig.pos <- sort(sample(n,length(x)))
        y[orig.pos] <- x
        warning("Vector to split is shorter than the number of desired bins! ",
            n-length(x)," 0s will be randomly distributed along the original ",
            "vector.",immediate.=TRUE)
        x <- y
    }
    bin.size <- floor(length(x)/n)
    dif <- length(x) - bin.size*n 
    bin.fac <- rep(bin.size,n)
    # Random bin increase size to avoid problems
    set.seed(seed)
    add <- sample(1:n,dif)
    bin.fac[add] <- bin.fac[add]+1
    f <- factor(rep(1:n,bin.fac))
    return(split(x,f))
}

readConfig <- function(input) {
    if (missing(input) || !file.exists(input))
        stop("File to read sample info from should be a valid existing text ",
            "file!")
    tab <- read.delim(input)
    samples <- as.character(tab[,1])
    if (length(samples) != length(unique(samples)))
        stop("Sample identifiers must be unique for each sample!")
    files <- as.character(tab[,3])
    formats <- as.character(tab[,4])
    if (any(!file.exists(files))) {
        bi <- which(!file.exists(files))
        stop("Input file ",files[bi]," does not exist! Please check paths...")
    }
    if (!all(formats %in% c("bam","bed")))
        stop("Input formats must be one of \"bam\", \"bed\"")
    output <- vector("list",nrow(tab))
    for (i in 1:nrow(tab)) {
        output[[i]]$id <- samples[i]
        output[[i]]$name <- as.character(tab[i,2])
        output[[i]]$file <- files[i]
        output[[i]]$format <- formats[i]
        output[[i]]$ranges <- NULL
        output[[i]]$coverage <- NULL
        output[[i]]$color <- NULL
    }
    return(output)
}

preprocessRanges <- function(input,preprocessParams,bamParams=NULL,rc=NULL) {
    hasRanges <- sapply(input,function(x) is.null(x$ranges))
    if (!any(hasRanges))
        return(input)
    switch(preprocessParams$normalize,
        none = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,params=p))
            },preprocessParams,bamParams,rc=rc)
            names(ranges) <- names(input)
            for (n in names(input))
                input[[n]]$ranges <- ranges[[n]]
        },
        linear = { # Same as none but will process after coverage calculation
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,params=p))
            },preprocessParams,bamParams,rc=rc)
            names(ranges) <- names(input)
            for (n in names(input))
                input[[n]]$ranges <- ranges[[n]]
        },
        downsample = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,params=p))
            },preprocessParams,bamParams,rc=rc)
            libSizes <- lengths(ranges)
            downto = min(libSizes)
            set.seed(preprocessParams$seed)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sort(sample(x,s)))
            },downto)
            names(downsampleIndex) <- names(input)
            for (n in names(input))
                input[[n]]$ranges <- ranges[[n]][downsampleIndex[[n]]]
        },
        sampleto = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,params=p))
            },preprocessParams,bamParams,rc=rc)
            set.seed(preprocessParams$seed)
            libSizes <- lengths(ranges)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sample(x,s))
            },preprocessParams$sampleTo)
            names(downsampleIndex) <- names(input)
            for (n in names(input)) {
                input[[n]]$ranges <- ranges[[n]][downsampleIndex[[n]]]
            }
        }
    )
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

readRanges <- function(input,format,sa,sq,params=NULL) {
    if (format=="bam")
        return(readBam(input,sa,sq,params))
    else if (format=="bed")
        return(readBed(input))
}

readBam <- function(bam,sa=c("keep","remove","split"),sq=0.75,params=NULL) {
    sa = sa[1]
    checkTextArgs("sa",sa,c("keep","remove","split"),multiarg=FALSE)
    checkNumArgs("sq",sq,"numeric",c(0,1),"botheq")
    switch(sa,
        keep = {
            return(trim(as(readGAlignments(file=bam),"GRanges")))
        },
        split = {
            return(trim(unlist(grglist(readGAlignments(file=bam)))))
        },
        remove = {
            reads <- trim(as(readGAlignments(file=bam),"GRanges"))
            qu <- quantile(width(reads),sq)
            rem <- which(width(reads)>qu)
            if (length(rem)>0)
                reads <- reads[-rem]
            message("  Excluded ",length(rem)," reads")
            return(reads)
        }
    )
}

readBed <- function(bed) {
    bed <- import.bed(sample.files[n],trackLine=FALSE,asRangedData=FALSE)
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!require(parallel)) || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.na(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            options(cores=ncores)
        }
    }
    if (m)
        return(mclapply(...,mc.cores=getOption("cores"),mc.set.seed=FALSE))
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
        binParams = {
            return(list(
                flankBinSize=0,
                regionBinSize=0,
                sumStat=c("mean","median"),
                smooth=TRUE,
                forceHeatmapBinning=TRUE,
                forcedBinSize=c(50,200)
            ))
        },
        preprocessParams = {
            return(list(
                normalize=c("none","linear","downsample","sampleto"),
                sampleTo=1e+6,
                spliceAction=c("keep","remove","split"),
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
                profile=TRUE,
                heatmap=TRUE,
                signalScale=c("natural","log2"),
                heatmapScale=c("each","common"),
                device=c("x11","png","jpg","tiff","bmp","pdf","ps"),
                outputDir=".",
                outputBase=NULL
            ))
        },
        saveParams = {
            return(list(
                ranges=TRUE,
                coverage=TRUE,
                profile=FALSE
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

getBiotypes <- function(org) {
    if (!(org %in% c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7",
        "pantro4","susscr3","tair10")))
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
