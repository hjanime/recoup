library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(circlize)

recoverProfile <- function(recoverObj,rc=NULL) {
    # Retrieve data
    input <- recoverObj$data
    design <- recoverObj$design
    opts <- recoverObj$plotopts
    
    # Create the x-axis breaks and labels
    width <- ncol(input[[1]]$profile)
    lb <- makeHorizontalAnnotation(width,opts)
    breaks <- lb$breaks
    labels <- lb$labels
    
    profileColors <- unlist(sapply(input,function(x) return(x$color)))
    if (!is.null(profileColors))
        names(profileColors) <- unlist(sapply(input,function(x) return(x$name)))
        
    if (is.null(design)) {
        # Create ggplot data
        profiles <- calcPlotProfiles(input,opts,rc)
        index <- 1:length(profiles[[1]][[1]])
        names <- sapply(input,function(x) return(x$name))
        position <- rep(index,length(input))
        signal <- unlist(lapply(profiles,function(x) return(x$profile)))
        ymin <- unlist(lapply(profiles,function(x) return(x$lower)))
        ymax <- unlist(lapply(profiles,function(x) return(x$upper)))
        condition <- rep(names,each=length(index))
        
        # Create ggplot data frame
        ggplot.data <- data.frame(
            Position=position,
            Signal=signal,
            Condition=factor(condition,levels=unique(condition)),
            ymin=ymin,
            ymax=ymax
        )
        
        # Create ggplot plot
        ggplot.plot <-
            ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,
                colour=Condition)) + 
            geom_line(size=1) +
            geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,colour=Condition,
                fill=Condition),alpha=0.3,size=0) +
            theme_bw() +
            xlab("\nPosition in bp") +
            ylab("Normalized average signal\n") +
            theme(title=element_text(size=14),
                axis.title.x=element_text(size=12,face="bold"),
                axis.title.y=element_text(size=12,face="bold"),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=12,face="bold"),
                legend.position="bottom") +
                scale_x_continuous(breaks=breaks,labels=labels)
                
        if (!is.null(profileColors))
            ggplot.plot <- ggplot.plot + 
                scale_fill_manual(values=profileColors) +
                scale_color_manual(values=profileColors)
    }
    else {
        message("Using provided design to facet the coverage profiles")
        subcovmat <- lapply(input,function(x,d) {
            splitter <- split(rownames(x$profile),as.list(d),sep="|",drop=TRUE)
            out <- vector("list",length(splitter))
            names(out) <- names(splitter)
            for (n in names(splitter))
                out[[n]] <- x$profile[splitter[[n]],,drop=FALSE]
            return(out)
        },design)
        subProfiles <- lapply(subcovmat,function(x,opts,rc) {
            return(calcDesignPlotProfiles(x,opts,rc))
        },opts,rc)
        
        designProfiles <- lapply(subProfiles,function(x) {
            d <- names(x) # Should be the "|" separated factor names
            dsplit <- strsplit(d,split="|",fixed=TRUE)
            # Replication factors
            pop <- sapply(x,function(x) return(length(x$profile)))
            tmp <- vector("list",length(x))
            for (i in 1:length(dsplit)) {
                if (length(dsplit[[i]])>1)
                    tmp[[i]] <- t(replicate(pop[i],dsplit[[i]]))
                else
                    tmp[[i]] <- as.matrix(rep(dsplit[[i]],pop[i]))
            }
            o <- list()
            o$profile <- unlist(lapply(x,function(y) {
                return(y$profile)
            }),use.names=FALSE)
            o$upper <- unlist(lapply(x,function(y) {
                return(y$upper)
            }),use.names=FALSE)
            o$lower <- unlist(lapply(x,function(y) {
                return(y$lower)
            }),use.names=FALSE)
            o$design <- do.call("rbind",tmp)
            return(o)
        })
        
        index <- 1:ncol(input[[1]]$profile)
        names <- sapply(input,function(x) return(x$name))
        faceter <- do.call("rbind",lapply(designProfiles,function(x) 
                return(x$design)))
        m <- length(which(!duplicated(faceter)))
        position <- rep(index,length(input)*m)
        condition <- rep(names,each=length(index)*m)
        signal <- unlist(lapply(designProfiles,function(x)
            return(x$profile)),use.names=FALSE)
        ymin <- unlist(lapply(designProfiles,function(x) 
            return(x$lower)),use.names=FALSE)
        ymax <- unlist(lapply(designProfiles,function(x) 
            return(x$upper)),use.names=FALSE)
        
        if (length(input)>1) { # Case where two factors max
            ggplot.data <- data.frame(
                Position=position,
                Signal=signal,
                Condition=factor(condition,levels=unique(condition)),
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1)
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            if (ncol(design)==2) {
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            
            ggplot.plot <-
                ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,
                    colour=Condition)) + 
                geom_line(size=1) +
                geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,
                    colour=Condition,fill=Condition),alpha=0.3,size=0) +
                theme_bw() +
                xlab("\nPosition in bp") +
                ylab("Normalized average signal\n") +
                theme(title=element_text(size=14),
                    axis.title.x=element_text(size=12,face="bold"),
                    axis.title.y=element_text(size=12,face="bold"),
                    axis.text.x=element_text(size=10,face="bold"),
                    axis.text.y=element_text(size=12,face="bold"),
                    strip.text.x=element_text(size=12,face="bold"),
                    strip.text.y=element_text(size=11,face="bold"),
                    legend.position="bottom",
                    panel.margin=unit(1,"lines"))  +
                scale_x_continuous(breaks=breaks,labels=labels)

            if (!is.null(profileColors))
                ggplot.plot <- ggplot.plot + 
                scale_fill_manual(values=profileColors) +
                scale_color_manual(values=profileColors)
            
            if (ncol(design)==1)
                ggplot.plot <- ggplot.plot + facet_grid(fac1~.)
            if (ncol(design)==2)
                ggplot.plot <- ggplot.plot + facet_grid(fac1~fac2)
        }
        else {
            ggplot.data <- data.frame(
                Position=position,
                Signal=signal,
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            }
            if (ncol(design)==2) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            if (ncol(design)==3) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
                ggplot.data$fac3 <- factor(as.character(faceter[,3]),
                    levels=unique(as.character(faceter[,3])))
            }
            
            ggplot.plot <-
                ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,
                    colour=Design)) + 
                geom_line(size=1) +
                geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,
                    colour=Design,fill=Design),alpha=0.3,size=0) +
                theme_bw() +
                xlab("\nPosition in bp") +
                ylab("Normalized average signal\n") +
                theme(title=element_text(size=14),
                    axis.title.x=element_text(size=12,face="bold"),
                    axis.title.y=element_text(size=12,face="bold"),
                    axis.text.x=element_text(size=10,face="bold"),
                    axis.text.y=element_text(size=12,face="bold"),
                    strip.text.x=element_text(size=12,face="bold"),
                    strip.text.y=element_text(size=11,face="bold"),
                    legend.position="bottom",
                    panel.margin=unit(1,"lines"))  +
                scale_x_discrete(breaks=breaks,labels=labels)

            if (ncol(design)==2)
                ggplot.plot <- ggplot.plot + facet_grid(fac2~.)
            if (ncol(design)==3)
                ggplot.plot <- ggplot.plot + facet_grid(fac2~fac3)
        }
   }
   return(ggplot.plot)
}

recoverHeatmap <- function(recoverObj,rc=NULL) {
    input <- recoverObj$data
    design <- recoverObj$design
    opts <- recoverObj$plotopts
    
    width <- ncol(input[[1]]$profile)
    labCol <- rep("",width)
    lb <- makeHorizontalAnnotation(width,opts)
    labCol[round(lb$breaks)] <- lb$labels
    
    if (opts$signalScale=="log2") {
        for (n in names(input)) {
            input[[n]]$profile <- input[[n]]$profile + 1
            input[[n]]$profile <- log2(input[[n]]$profile)
        }
    }
    
    haCol <- HeatmapAnnotation(cn=function(index) {
        width <- ncol(input[[1]]$profile)
        labCol <- rep("",width)
        lb <- makeHorizontalAnnotation(width,opts)
        labCol[round(lb$breaks)] <- lb$labels
        grid.text(labCol,(1:width)/width,1,vjust=1,
            gp=gpar(cex=0.8))
    })
    
    # Here we need to apply the orderBy directives. If orderBy$what is 
    # something else that hierarchical clustering, parameter control in main
    # has already fixed incompatibilities by turning off clustering. If it
    # is clustering, it's turned on, so we need to decide which is the main
    # profile.
    # The same ordering ideas must be applied to the design elements... We must
    # reorder the rownames of the design matrix in such a way as as to be ready 
    # when fed to Heatmap function.
    if (is.null(design))
        sorter <- orderProfiles(input,opts,rc=rc)
    else
        sorter <- orderProfilesByDesign(input,design,opts,rc=rc)
    
    colorFuns <- vector("list",length(input))
    names(colorFuns) <- names(input)
    profileColors <- unlist(sapply(input,function(x) return(x$color)))
    if (opts$heatmapScale=="each") {
        for (n in names(colorFuns)) {
            qs <- c(0.95,0.96,0.97,0.98,0.99,0.995,0.999)
            pos <- 1
            sup <- quantile(input[[n]]$profile,qs[pos])
            while(sup==0) {
                pos <- pos + 1
                sup <- quantile(input[[n]]$profile,qs[pos])
                if (sup!=0)
                    break
            }
            if (!is.null(profileColors))
                colorFuns[[n]] <- colorRamp2(c(0,sup),c("white",
                    input[[n]]$color))
            else
                colorFuns[[n]] <- colorRamp2(c(0,sup),c("white","red2"))
        }
    }
    else if (opts$heatmapScale=="common") {
        sups <- unlist(sapply(input,function(x) {
            return(quantile(x$profile,0.95))
        }))
        sup <- max(sups)
        for (n in names(colorFuns)) {
            if (!is.null(profileColors))
                colorFuns[[n]] <- colorRamp2(c(0,sup), c("white",
                    input[[n]]$color))
            else
                colorFuns[[n]] <- colorRamp2(c(0,sup), c("white","red2"))
        }
    }
    
    if (is.null(design)) {
        hParams <- opts$complexHeatmapParams$main
        hmList <- NULL
        for (n in names(input)) {
            colnames(input[[n]]$profile) <- labCol
            hmList <- hmList + 
                Heatmap(
                    input[[n]]$profile[sorter$ix,],
                    name=input[[n]]$name,
                    cluster_rows=hParams$cluster_rows,
                    cluster_columns=hParams$cluster_columns,
                    column_title_gp=hParams$column_title_gp,
                    show_row_names=hParams$show_row_names,
                    show_column_names=hParams$show_column_names,
                    col=colorFuns[[n]],
                    column_title=paste(input[[n]]$name,"signal"),
                    heatmap_legend_param=hParams$heatmap_legend_param,
                    bottom_annotation=haCol
                )
        }
    }
    else {
        message("Using provided design to facet the coverage profiles")
        hParams <- opts$complexHeatmapParams$group
        hmList <- NULL
        for (n in names(input)) {
            colnames(input[[n]]$profile) <- labCol
            hmList <- hmList + 
                Heatmap(
                    input[[n]]$profile,
                    name=input[[n]]$name,
                    cluster_rows=hParams$cluster_rows,
                    cluster_columns=hParams$cluster_columns,
                    column_title_gp=hParams$column_title_gp,
                    show_row_names=hParams$show_row_names,
                    show_column_names=hParams$show_column_names,
                    col=colorFuns[[n]],
                    column_title=paste(input[[n]]$name,"signal"),
                    heatmap_legend_param=hParams$heatmap_legend_param,
                    bottom_annotation=haCol,
                    split=design,
                    row_order=sorter,
                    row_title_gp=hParams$row_title_gp,
                    gap=hParams$gap
                )
       }
   }
   return(hmList)
}

calcPlotProfiles <- function(input,opts,rc) {
    if (opts$binParams$smooth)
        profiles <- cmclapply(input,function(x,avgfun,scale) {
            if (scale=="log2") {
                x$profile <- x$profile + 1
                x$profile <- log2(x$profile)
            }
            o <- list()
            fit <- smooth.spline(apply(x$profile,2,avgfun))
            ci <- ssCI(fit)
            o$profile <- fit$y
            o$upper <- ci$upper
            o$lower <- ci$lower
            return(o)
        },opts$binParams$sumStat,opts$signalScale,rc=rc)
    else
        profiles <- cmclapply(input,function(x,avgfun,scale) {
            if (scale=="log2") {
                x$profile <- x$profile + 1
                x$profile <- log2(x$profile)
            }
            o <- list()
            o$profile <- apply(x$profile,2,avgfun)
            varfun <- ifelse(avgfun=="mean","sd","mad")
            va <- apply(x,2,varfun)
            o$upper <- o$profile + va
            o$lower <- o$profile - va
            return(o)
        },opts$binParams$sumStat,opts$signalScale,rc=rc)
    return(profiles)
}

calcDesignPlotProfiles <- function(covmat,opts,rc) {
    if (opts$binParams$smooth)
        profiles <- cmclapply(covmat,function(x,avgfun,scale) {
            if (scale=="log2") {
                x <- x + 1
                x <- log2(x)
            }
            o <- list()
            fit <- smooth.spline(apply(x,2,avgfun))
            ci <- ssCI(fit)
            o$profile <- fit$y
            o$upper <- ci$upper
            o$lower <- ci$lower
            return(o)
        },opts$binParams$sumStat,opts$signalScale,rc=rc)
    else
        profiles <- cmclapply(covmat,function(x,avgfun,scale) {
            if (scale=="log2") {
                x <- x + 1
                x <- log2(x)
            }
            o <- list()
            o$profile <- apply(x,2,avgfun)
            varfun <- ifelse(avgfun=="mean","sd","mad")
            va <- apply(x,2,varfun)
            o$upper <- o$profile + va
            o$lower <- o$profile - va
            return(o)
        },opts$binParams$sumStat,opts$signalScale,rc=rc)
    return(profiles)
}

orderProfiles <- function(input,opts,rc=NULL) {
    refh <- 1
    if (length(grep("^(sum|max)",opts$orderBy$what,perl=TRUE))>0) {
        nc <- nchar(opts$orderBy$what)
        rh <- substr(opts$orderBy$what,nc,nc)
        if (rh!="a") {
            rh <- suppressWarnings(as.numeric(rh))
            if (is.na(rh))
                warning("Reference profile for heatmap ordering not ",
                    "recognized! Using the 1st...",immediate.=TRUE)
            else
                refh <- rh
        }
        else
            refh <- 0 # Flag to indicate order by sum/max of all profiles
    }
    byMax <- bySum <- FALSE
    sorter <- list(ix=1:nrow(input[[1]]$profile))
    if (length(grep("^sum",opts$orderBy$what,perl=TRUE))>0)
        bySum <- TRUE
    if (length(grep("^max",opts$orderBy$what,perl=TRUE))>0)
        byMax <- TRUE
    
    if (bySum) {
        if (refh==0) {
            #theBigMatrix <- do.call("cbind",lapply(input,function(x) {
            #   return(x$profile)
            #}))
            #theSum <- apply(theBigMatrix,1,sum)
            tmp <- do.call("cbind",lapply(input,function(x,rc) {
                s <- cmclapply(x$coverage,function(y) {
                    if (is.list(y))
                        return(sum(y$center))
                    else
                        return(sum(y))
                },rc=rc)
                return(unlist(s))
            },rc=rc))
            theVal <- apply(tmp,1,sum)
            names(theVal) <- rownames(input[[1]]$profile)
        }
        else {
            theVal <- unlist(cmclapply(input[[refh]]$coverage,function(y) {
                if (is.list(y))
                    return(sum(y$center))
                else
                    return(sum(y))
            },rc=rc))
            names(theVal) <- rownames(input[[refh]]$profile)
        }
        if (opts$orderBy$order=="descending")
            sorter <- sort(theVal,decreasing=TRUE,index.return=TRUE)
        else
            sorter <- sort(theVal,index.return=TRUE)
    }
    if (byMax) {
        if (refh==0) {
            tmp <- do.call("cbind",lapply(input,function(x,rc) {
                s <- cmclapply(x$coverage,function(y) {
                    if (is.list(y))
                        y <- y$center
                    m <- max(y)
                    mp <- which(y==m)
                    if (length(mp)>1)
                        return(as.numeric(y[sample(mp,1)]))
                    else
                        return(as.numeric(y[mp]))
                },rc=rc)
                return(unlist(s))
            },rc=rc))
            theVal <- apply(tmp,1,max)
            names(theVal) <- rownames(input[[1]]$profile)
        }
        else {
            theVal <- unlist(cmclapply(input[[refh]]$coverage,function(y) {
                if (is.list(y))
                    y <- y$center
                m <- max(y)
                mp <- which(y==m)
                if (length(mp)>1)
                    return(as.numeric(y[sample(mp,1)]))
                else
                    return(as.numeric(y[mp]))
            },rc=rc))
            names(theVal) <- rownames(input[[refh]]$profile)
        }
        if (opts$orderBy$order=="descending")
            sorter <- sort(theVal,decreasing=TRUE,index.return=TRUE)
        else
            sorter <- sort(theVal,index.return=TRUE)
    }
    return(sorter)
}

orderProfilesByDesign <- function(input,design,opts,rc=NULL) {
    splitter <- split(1:nrow(design),design,drop=TRUE)
    sortlist <- sorter <- vector("list",length(splitter))
    names(sortlist) <- names(sorter) <- names(splitter)
    
    for (n in names(splitter)) {
        S <- splitter[[n]]
        subcov <- lapply(input,function(x,s) {
            if (!is.null(x$coverage$center))
                return(x$coverage$center[s])
            else
                return(x$coverage[s])
        },S)
            
        refh <- 1
        if (length(grep("^(sum|max)",opts$orderBy$what,perl=TRUE))>0) {
            nc <- nchar(opts$orderBy$what)
            rh <- substr(opts$orderBy$what,nc,nc)
            if (rh!="a") {
                rh <- suppressWarnings(as.numeric(rh))
                if (is.na(rh))
                    warning("Reference profile for heatmap ordering not ",
                        "recognized! Using the 1st...",immediate.=TRUE)
                else
                    refh <- rh
            }
            else
                refh <- 0 # Flag to indicate order by sum/max of all profiles
        }
        byMax <- bySum <- FALSE
        sortlist[[n]] <- 1:length(S)
        if (length(grep("^sum",opts$orderBy$what,perl=TRUE))>0)
            bySum <- TRUE
        if (length(grep("^max",opts$orderBy$what,perl=TRUE))>0)
            byMax <- TRUE
        
        if (bySum) {
            if (refh==0) {
                tmp <- do.call("cbind",lapply(subcov,function(x,rc) {
                    s <- cmclapply(x,function(y) {
                        return(sum(y))
                    },rc=rc)
                    return(unlist(s))
                },rc=rc))
                theVal <- apply(tmp,1,sum)
            }
            else {
                theVal <- unlist(cmclapply(subcov[[refh]],function(y) {
                    return(sum(y))
                },rc=rc))
            }
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(theVal,decreasing=TRUE,index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(theVal,index.return=TRUE)$ix
        }
        if (byMax) {
            if (refh==0) {
                tmp <- do.call("cbind",lapply(subcov,function(x,rc) {
                    s <- cmclapply(x,function(y) {
                        m <- max(y)
                        mp <- which(y==m)
                        if (length(mp)>1)
                            return(as.numeric(y[sample(mp,1)]))
                        else
                            return(as.numeric(y[mp]))
                    },rc=rc)
                    return(unlist(s))
                },rc=rc))
                theVal <- apply(tmp,1,max)
            }
            else {
                theVal <- unlist(cmclapply(subcov[[refh]],function(y) {
                    m <- max(y)
                    mp <- which(y==m)
                    if (length(mp)>1)
                        return(as.numeric(y[sample(mp,1)]))
                    else
                        return(as.numeric(y[mp]))
                },rc=rc))
            }
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(theVal,decreasing=TRUE,index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(theVal,index.return=TRUE)$ix
        }
        sorter[[n]] <- S[sortlist[[n]]]
   }
   return(unlist(sorter))
}

makeHorizontalAnnotation <- function(width,opts) {
    fl <- opts$flank
    fl[1] <- -fl[1]
    edgeLabels <- paste(round(fl/1000,1),"kb",sep="")
    switch(opts$region,
        tss = {
            midLabels <- "TSS"
            #breaks <- round(c(width/8,width/2,width-width/8),1)
            breaks <- c(1,round(width/2,1),width)
        },
        tes = {
            midLabels <- "TES"
            #breaks <- round(c(width/8,width/2,width-width/8),1)
            breaks <- c(1,round(width/2,1),width)
        },
        genebody = {
            midLabels <- c("TSS","TES")
            if (opts$binParams$flankBinSize==0)
                breaks <- c(
                    1, #round(abs(opts$flank[1])/8,1),
                    abs(opts$flank[1]),
                    width-opts$flank[2],
                    width #round(width-opts$flank[2]/8)
                )
            else
                breaks <- c(
                    1, #round(opts$binParams$flankBinSize/8,1),
                    opts$binParams$flankBinSize,
                    width-opts$binParams$flankBinSize,
                    width #round(width-opts$binParams$flankBinSize/8)
                )
        },
        custom = {
            if (opts$customOne) {
                midLabels <- "Center"
                breaks <- round(c(width/8,width/2,width-width/8),1)
            }
            else {
                midLabels <- c("Start","End")
                if (opts$binParams$flankBinSize==0)
                    breaks <- c(
                        1, #round(abs(opts$flank[1])/8,1),
                        abs(opts$flank[1]),
                        width-opts$flank[2],
                        width #round(width-opts$flank[2]/8)
                    )
                else
                    breaks <- c(
                        1, #round(opts$binParams$flankBinSize/8,1),
                        opts$binParams$flankBinSize,
                        width-opts$binParams$flankBinSize,
                        width #round(width-opts$binParams$flankBinSize/8)
                    )
            }
        }
    )
    labels <- c(edgeLabels[1],midLabels,edgeLabels[2])
    return(list(breaks=breaks,labels=labels))
}

graphicsOpen <- function(o,f,...) {
    if(o!="x11" && is.null(f))
        stop("Please specify an output file name for your plot")
    switch(o,
        x11 = { dev.new(...) },
        pdf = { pdf(file=f,pointsize=10,...) },
        ps = { postscript(file=f,pointsize=10,...) },
        png = { png(filename=f,pointsize=12,...) },
        jpg = { jpeg(filename=f,pointsize=12,quality=100,...) },
        bmp = { bmp(filename=f,pointsize=12,...) },
        tiff = { tiff(filename=f,pointsize=12,...) }
    )
}

graphicsClose <- function(o) {
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    if (o!="x11")
        dev.off()
}
