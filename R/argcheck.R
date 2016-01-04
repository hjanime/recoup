checkTextArgs <- function(arg.name,arg.value,arg.list,multiarg=FALSE) {
    if (multiarg) {
        arg.value <- tolower(arg.value)
        if (!all(arg.value %in% arg.list))
            stop("\"",arg.name,"\""," parameter must be one or more of ",
                paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
    }
    else {
        arg.save <- arg.value[1]
        arg.value <- tolower(arg.value[1])
        # An exception must be added for annotation because it can be an external 
        # file too
        if (!(arg.value %in% arg.list))
            stop("\"",arg.name,"\""," parameter must be one of ",
                paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
    }
}

checkNumArgs <- function(arg.name,arg.value,arg.type,arg.bounds,direction) {
    switch(arg.type,
        numeric = {
            if (!is.numeric(arg.value))
                stop("\"",arg.name,"\"",
                    " parameter must be a numeric value!")
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric ","value larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric value larger than ",arg.bounds[1],
                                " and smaller than ",arg.bounds[2],"!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric value greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric value lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric value greater than or equal to ",
                                arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be a ",
                                "numeric value lower than or equal to ",
                                arg.bounds[1],"!")
                    }
                )
            }
        },
        integer = {
            if (!is.integer(arg.value))
                stop("\"",arg.name,"\""," parameter must be an integer!")
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer greater than or equal to ",
                                arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stop("\"",arg.name,"\""," parameter must be an ",
                                "integer lower than or equal to ",
                                arg.bounds[1],"!")
                    }
                )
            }
        }
    )
}

checkFileArgs <- function(arg.name,arg.value) {
    if (!file.exists(arg.value))
        stop("\"",arg.name,"\""," parameter must be an existing file!")
}

validateListArgs <- function(what,arg.list) {
    if (!(what %in% c("binParams","preprocessParams","selector",
        "strandedParams","saveParams","plotParams","kmParams")))
        stop("Input list type for validation must be one of \"binParams\", ",
            "\"preprocessParams\", \"selector\", \"strandedParams\", ",
            "\"saveParams\", \"plotParams\", \"kmParams\"")
    if (!is.list(arg.list))
        stop(what," argument must be a list!")
    switch(what,
        binParams = {
            valid.1 <- names(arg.list) %in% c("flankBinSize","regionBinSize",
                "sumStat","smooth","forceHeatmapBinning","forcedBinSize")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        },
        preprocessParams = {
            valid.1 <- names(arg.list) %in% c("normalize","sampleTo",
                "spliceAction","seed")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            #if (length(arg.list)>0) {
            #    for (n in names(arg.list)) {
            #        switch(n,
            #            normalize = {
            #                valid.2 <- names(arg.list[[n]]) %in% c("none",
            #                    "linear","downsample","sampleto")
            #                not.valid.2 <- which(!valid.2)
            #            },
            #            spliceAction = {
            #                valid.2 <- names(arg.list[[n]]) %in%
            #                    c("keep","remove","split")
            #                not.valid.2 <- which(!valid.2)
            #            }
            #        )
            #        if (length(not.valid.2)>0) {
            #            warning("The following ",what," sub-argument names ",
            #                "are invalid and will be ignored: ",
            #                paste(names(arg.list[[n]])[not.valid.2],
            #                    collapse=", "))
            #            arg.list[[n]][not.valid.2] <- NULL
            #        }
            #    }
            #}
        },
        selector = {
            valid.1 <- names(arg.list) %in% c("id","biotype","exonType")
            not.valid.1 <- which(!valid)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        },
        strandedParams = {
            valid.1 <- names(arg.list) %in% c("strand","ignoreStrand")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        },
        saveParams = {
            valid.1 <- names(arg.list) %in% c("ranges","coverage","profile")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        },
        plotParams = {
            valid.1 <- names(arg.list) %in% c("profile","heatmap","device",
                "heatmapScale","outputDir","outputBase")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        },
        kmParams = {
            valid.1 <- names(arg.list) %in% c("k","nstart","algorithm",
                "reference","iterMax","seed")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
        }
    )
    return(arg.list)
}

checkInput <- function(input) {
    # Input must have  id, file and format
    for (i in 1:length(input)) {
        if (is.null(input[[i]]$id))
            stop("All input list members must have an id field! Member ",i,
                " does not have one.")
        if (is.null(input[[i]]$file))
            stop("All input list members must have a file field! Member ",i,
                " does not have one.")
        if (is.null(input[[i]]$format))
            stop("All input list members must have a format field! Member ",i,
                " does not have one.")
    }
}
