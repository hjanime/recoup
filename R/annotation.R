require(biomaRt)

buildAnnotationStore <- function(organisms,sources,
    home=file.path(path.expand("~"),".recover"),rc=NULL) {
    # Load required packages
    require(GenomicRanges)

    if (missing(organisms))
        organisms <- c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7",
            "pantro4","susscr3")
    if (missing(sources))
        sources <- c("ensembl","ucsc","refseq")
    
    checkTextArgs("organisms",organisms,c("hg18","hg19","hg38","mm9","mm10",
        "rn5","dm3","danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("sources",sources,c("ensembl","ucsc","refseq"),multiarg=FALSE)
    
    # Retrieve gene annotations
    for (s in sources) {
        for (o in organisms) {
            message("Retrieving gene annotation for ",o," from ",s)
            ann <- getAnnotation(o,"gene",refdb=s,rc=rc)
            gene <- makeGRangesFromDataFrame(
                df=ann,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome"
            )
            store.path <- file.path(home,s,o)
            if (!dir.exists(store.path))
                dir.create(store.path,recursive=TRUE,mode="0755")
            save(gene,file=file.path(store.path,"gene.rda"),compress=TRUE)
        }
    }

    # Retrieve exon annotations
    for (s in sources) {
        for (o in organisms) {
            message("Retrieving exon annotation for ",o," from ",s)
            ann <- getAnnotation(o,"exon",refdb=s,rc=rc)
            ann.gr <- makeGRangesFromDataFrame(
                df=ann,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome"
            )
            exon <- split(ann.gr,ann.gr$gene_id)
            store.path <- file.path(home,s,o)
            if (!dir.exists(store.path))
                dir.create(store.path,recursive=TRUE,mode="0755")
            save(exon,file=file.path(store.path,"exon.rda"),compress=TRUE)
            # Then summarize the exons and write again with type sum_exon
            message("Merging exons for ",o," from ",s)
            ann.gr <- reduceExons(ann.gr,rc=rc)
            names(ann.gr) <- as.character(ann.gr$exon_id)
            sexon <- split(ann.gr,ann.gr$gene_id)
            store.path <- file.path(home,s,o)
            if (!dir.exists(store.path))
                dir.create(store.path,recursive=TRUE,mode="0755")
            save(sexon,file=file.path(store.path,"summarized_exon.rda"),
                compress=TRUE)
        }
    }
}

#' Merges exons to create a unique set of exons for each gene
#'
#' This function uses the \code{"reduce"} function of IRanges to construct virtual
#' unique exons for each gene, so as to avoid inflating the read counts for each
#' gene because of multiple possible transcripts. If the user wants transcripts
#' instead of genes, they should be supplied to the original annotation table.
#'
#' @param gr a GRanges object created from the supplied annotation (see also the
#' \code{\link{read2count}} and \code{\link{get.annotation}} functions.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not.
#' @return A GRanges object with virtual merged exons for each gene/transcript.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(GenomicRanges)
#' ann <- get.annotation("mm9","exon")
#' gr <- makeGRangesFromDataFrame(
#'  df=ann,
#'  keep.extra.columns=TRUE,
#'  seqnames.field="chromosome"
#' )
#' re <- reduce.exons(gr)
#'}
reduceExons <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    red.list <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            exon_id=paste(x,"MEX",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    return(do.call("c",red.list))
}

#' Annotation downloader
#'
#' This function connects to the EBI's Biomart service using the package biomaRt
#' and downloads annotation elements (gene co-ordinates, exon co-ordinates, gene
#' identifications, biotypes etc.) for each of the supported organisms. See the
#' help page of \code{\link{metaseqr}} for a list of supported organisms. The
#' function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb the online source to use to fetch annotation. It can be
#' \code{"ensembl"} (default), \code{"ucsc"} or \code{"refseq"}. In the later two
#' cases, an SQL connection is opened with the UCSC public databases.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not. It is used in the case of \code{type="exon"} to process
#' the return value of the query to the UCSC Genome Browser database.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to Ensembl gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.annotation("hg19","gene","ensembl")
#' mm9.exons <- get.annotation("mm9","exon","ucsc")
#'}
getAnnotation <- function(org,type,refdb="ensembl",rc=NULL) {
    org <- tolower(org)
    switch(refdb,
        ensembl = { return(getEnsemblAnnotation(org,type)) },
        ucsc = { return(getUcscAnnotation(org,type,refdb,rc=rc)) },
        refseq = { return(getUcscAnnotation(org,type,refdb,rc=rc)) }
    )
}

#' Ensembl annotation downloader
#'
#' This function connects to the EBI's Biomart service using the package biomaRt
#' and downloads annotation elements (gene co-ordinates, exon co-ordinates, gene
#' identifications, biotypes etc.) for each of the supported organisms. See the
#' help page of \code{\link{metaseqr}} for a list of supported organisms. The
#' function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to Ensembl gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.ensembl.annotation("hg19","gene")
#' mm9.exons <- get.ensembl.annotation("mm9","exon")
#'}
getEnsemblAnnotation <- function(org,type) {
    dat <- "ENSEMBL_MART_ENSEMBL"
    mart <- tryCatch({
            useMart(biomart=dat,host=getHost(org),dataset=getDataset(org))
        },
        error=function(e) {
            useMart(biomart=dat,host=getAltHost(org),
            dataset=getDataset(org))
        },
        finally={})
    chrs.exp <- paste(getValidChrs(org),collapse="|")
    if (type=="gene") {
        bm <- getBM(attributes=getGeneAttributes(org),mart=mart)
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$start_position,
            end=bm$end_position,
            gene_id=bm$ensembl_gene_id,
            gc_content=bm$percentage_gc_content,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","mm9","tair10")) bm$external_gene_id 
                else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        bm <- getBM(attributes=getExonAttributes(org),mart=mart)
        if (org == "hg19") {
            disp("  Bypassing problem with hg19 Ensembl combined gene-exon ",
                "annotation... Will take slightly longer...")
            bmg <- getBM(attributes=getGeneAttributes(org),mart=mart)
            gene_name <- bmg$external_gene_name
            names(gene_name) <- bmg$ensembl_gene_id
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=gene_name[bm$ensembl_gene_id],
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
        }
        else
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=if (org %in% c("hg18","mm9","tair10")) 
                    bm$external_gene_id else bm$external_gene_name,
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
    }
    ann <- ann[order(ann$chromosome,ann$start),]
    ann <- ann[grep(chrs.exp,ann$chromosome),]
    ann$chromosome <- as.character(ann$chromosome)
    return(ann)
}

#' UCSC/RefSeq annotation downloader
#'
#' This function connects to the UCSC Genome Browser public database and downloads
#' annotation elements (gene co-ordinates, exon co-ordinates, gene identifications
#' etc.) for each of the supported organisms, but using UCSC instead of Ensembl.
#' See the help page of \code{\link{metaseqr}} for a list of supported organisms.
#' The function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb either \code{"ucsc"} or \code{"refseq"}.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not. It is used in the case of \code{type="exon"} to process
#' the return value of the query to the UCSC Genome Browser database.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to UCSC or RefSeq gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M. Note also that as the UCSC databases do not contain
#' biotype classification like Ensembl, this will be returned as \code{NA} and
#' as a result, some quality control plots will not be available.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- getUcscAnnotation("hg19","gene","ucsc")
#' mm9.exons <- getUcscAnnotation("mm9","exon")
#'}
getUcscAnnotation <- function(org,type,refdb="ucsc",rc=NULL) {
    if (!require(RMySQL)) {
        rmysql.present <- FALSE
        warnwrap("R package RMySQL is not present! Annotation will be ",
            "retrieved by downloading temporary files from UCSC and the usage
            of a temporary SQLite database...")
    }
    else
        rmysql.present <- TRUE
    if (!require(RSQLite))
        stop("R package RSQLite is required to use annotation from UCSC!")

    valid.chrs <- getValidChrs(org)
    chrs.exp <- paste("^",paste(valid.chrs,collapse="$|^"),"$",sep="")

    db.org <- getUcscOrganism(org)
    if (rmysql.present) {
        db.creds <- getUcscCredentials()
        drv <- dbDriver("MySQL")
        con <- dbConnect(drv,user=db.creds[2],password=NULL,dbname=db.org,
            host=db.creds[1])
        query <- getUcscQuery(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    else {
        # This should return the same data frame as the db query
        tmp.sqlite <- getUcscDbl(org,type,refdb)
        drv <- dbDriver("SQLite")
        con <- dbConnect(drv,dbname=tmp.sqlite)
        query <- getUcscQuery(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    if (type=="gene") {
        ann <- raw.ann
        ann <- ann[grep(chrs.exp,ann$chromosome,perl=TRUE),]
        ann$chromosome <- as.character(ann$chromosome)
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        raw.ann <- raw.ann[grep(chrs.exp,raw.ann$chromosome,perl=TRUE),]
        ex.list <- cmclapply(as.list(1:nrow(raw.ann)),function(x,d,s) {
            r <- d[x,]
            starts <- as.numeric(strsplit(r[,"start"],",")[[1]])
            ends <- as.numeric(strsplit(r[,"end"],",")[[1]])
            nexons <- length(starts)
            ret <- data.frame(
                rep(r[,"chromosome"],nexons),
                starts,ends,
                paste(r[,"exon_id"],"_e",1:nexons,sep=""),
                rep(r[,"strand"],nexons),
                rep(r[,"gene_id"],nexons),
                rep(r[,"gene_name"],nexons),
                rep(r[,"biotype"],nexons)
            )
            names(ret) <- names(r)
            rownames(ret) <- ret$exon_id
            ret <- makeGRangesFromDataFrame(
                df=ret,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome",
                seqinfo=s
            )
            return(ret)
        },raw.ann,valid.chrs,rc=rc)
        tmp.ann <- do.call("c",ex.list)
        ann <- data.frame(
            chromosome=as.character(seqnames(tmp.ann)),
            start=start(tmp.ann),
            end=end(tmp.ann),
            exon_id=as.character(tmp.ann$exon_id),
            gene_id=as.character(tmp.ann$gene_id),
            strand=as.character(strand(tmp.ann)),
            gene_name=as.character(tmp.ann$gene_name),
            biotype=tmp.ann$biotype
        )
        rownames(ann) <- ann$exon_id
    }
    
    gc.content <- getGcContent(ann,org)
    ann$gc_content <- gc.content
    ann <- ann[order(ann$chromosome,ann$start),]
    return(ann)
}

#' Return a named vector of GC-content for each genomic region
#'
#' Returns a named numeric vector (names are the genomic region names, e.g. genes)
#' given a data frame which can be converted to a GRanges object (e.g. it has at
#' least chromosome, start, end fields). This function works best when the input
#' annotation data frame has been retrieved using one of the SQL queries generated
#' from \code{\link{getUcscQuery}}, used in \code{\link{getUcscAnnotation}}.
#'
#' @param ann a data frame which can be converted to a GRanges object, that means
#' it has at least the chromosome, start, end fields. Preferably, the output of
#' \code{link{getUcscAnnotation}}.
#' @param org one of metaseqR supported organisms.
#' @return A named numeric vector.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' ann <- getUcscAnnotation("mm9","gene","ucsc")
#' gc <- getGcContent(ann,"mm9")
#'}
getGcContent <- function(ann,org) {
    if (missing(ann))
        stopwrap("A valid annotation data frame must be provided in order to ",
            "retrieve GC-content.")
    org <- tolower(org[1])
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    # Convert annotation to GRanges
    disp("Converting annotation to GenomicRanges object...")
    if (packageVersion("GenomicRanges")<1.14)
        ann.gr <- GRanges(
            seqnames=Rle(ann[,1]),
            ranges=IRanges(start=ann[,2],end=ann[,3]),
            strand=Rle(ann[,6]),
            name=as.character(ann[,4])
        )
    else
        ann.gr <- makeGRangesFromDataFrame(
            df=ann,
            keep.extra.columns=TRUE,
            seqnames.field="chromosome"
        )
    bsg <- loadBsGenome(org)
    disp("Getting DNA sequences...")
    seqs <- getSeq(bsg,names=ann.gr)
    disp("Getting GC content...")
    freq.matrix <- alphabetFrequency(seqs,as.prob=TRUE,baseOnly=TRUE)
    gc.content <- apply(freq.matrix,1,function(x) round(100*sum(x[2:3]),
        digits=2))
    names(gc.content) <- as.character(ann[,4])
    return(gc.content)
}

#' Return a proper formatted organism alias
#'
#' Returns the proper UCSC Genome Browser database organism alias based on what is
#' given to metaseqR. Internal use.
#'
#' @return A proper organism alias.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' org <- getUcscOrganism("danrer7")
#'}
getUcscOrganism <- function(org) {
    switch(org,
        hg18 = { return("hg18") },
        hg19 = { return("hg19") },
        hg38 = { return("hg38") },
        mm9 = { return("mm9") },
        mm10 = { return("mm10") },
        rn5 = { return("rn5") },
        dm3 = { return("dm3") },
        danrer7 = { return("danRer7") },
        pantro4 = { return("panTro4") },
        susscr3 = { return("susScr3") }
    )
}

#' Return a proper formatted BSgenome organism name
#'
#' Returns a properly formatted BSgenome package name according to metaseqR's
#' supported organism. Internal use.
#'
#' @return A proper BSgenome package name.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' bs.name <- getBsOrganism("hg18")
#'}
getBsOrganism <- function(org) {
    switch(org,
        hg18 = {
            return("BSgenome.Hsapiens.UCSC.hg18")
        },
        hg19 = {
            return("BSgenome.Hsapiens.UCSC.hg19")
        },
        hg38 = {
            return("BSgenome.Hsapiens.UCSC.hg38")
        },
        mm9 = {
            return("BSgenome.Mmusculus.UCSC.mm9")
        },
        mm10 = {
            return("BSgenome.Mmusculus.UCSC.mm10")
        },
        rn5 = {
            return("BSgenome.Rnorvegicus.UCSC.rn5")
        },
        dm3 = {
            return("BSgenome.Dmelanogaster.UCSC.dm3")
        },
        danrer7 = {
            return("BSgenome.Drerio.UCSC.danRer7")
        },
        pantro4 = {
            stopwrap("panTro4 is not yet supported by BSgenome! Please use ",
                "Ensembl as annoation source.")
        },
        susscr3 = {
            return("BSgenome.Sscrofa.UCSC.susScr3")
        }
    )
}

#' Loads (or downloads) the required BSGenome package
#'
#' Retrieves the required BSgenome package when the annotation source is \code{"ucsc"}
#' or \code{"refseq"}. These packages are required in order to estimate the
#' GC-content of the retrieved genes from UCSC or RefSeq.
#'
#' @param org one of \code{\link{metaseqr}} supported organisms.
#' @return The BSgenome object for the requested organism.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' bs.obj <- loadBsGenome("mm9")
#'}
loadBsGenome <- function(org) {
    if (!require(BiocInstaller))
        stopwrap("The Bioconductor package BiocInstaller is required to ",
            "proceed!")
    if (!require(BSgenome))
        stopwrap("The Bioconductor package BSgenome is required to ",
            "proceed!")
    bs.org <- getBsOrganism(org)
    if (bs.org %in% installed.genomes())
        bs.obj <- getBSgenome(getUcscOrganism(org))
    else {
        biocLite(bs.org)
        bs.obj <- getBSgenome(getUcscOrganism(org))
    }
    return(bs.obj)
}

#' Biotype converter
#'
#' Returns biotypes as character vector. Internal use.
#'
#' @param a the annotation data frame (output of \code{\link{get.annotation}}).
#' @return A character vector of biotypes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg18.genes <- get.annotation("hg18","gene")
#' hg18.bt <- getBiotypes(hg18.genes)
#'}
getBiotypes <- function(a) {
    return(as.character(unique(a$biotype)))
}

#' Annotation downloader helper
#'
#' Returns the appropriate Ensembl host address to get different versions of
#' annotation from. Internal use.
#'
#' @param org the organism for which to return the host address.
#' @return A string with the host address.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' mm9.host <- getHost("mm9")
#'}
getHost <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("www.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("www.ensembl.org") },
        rn5 = { return("www.ensembl.org") },
        dm3 = { return("www.ensembl.org") },
        danrer7 = { return("www.ensembl.org") },
        pantro4 = { return("www.ensembl.org") },
        susscr3 = { return("www.ensembl.org") },
        tair10 = { return("www.biomart.org") },
        bmori2 = { return("metazoa.ensembl.org") }
    )
}

#' Annotation downloader helper
#'
#' Returns the appropriate Ensembl host address to get different versions of
#' annotation from (alternative hosts). Internal use.
#'
#' @param org the organism for which to return the host address.
#' @return A string with the host address.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' mm9.host <- getAltHost("mm9")
#'}
getAltHost <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("uswest.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("uswest.ensembl.org") },
        rn5 = { return("uswest.ensembl.org") },
        dm3 = { return("uswest.ensembl.org") },
        danrer7 = { return("uswest.ensembl.org") },
        pantro4 = { return("uswest.ensembl.org") },
        susscr3 = { return("uswest.ensembl.org") },
        tair10 = { return("www.biomart.org") },
        bmori2 = { return("metazoa.ensembl.org") }
    )
}

#' Annotation downloader helper
#'
#' Returns a dataset (gene or exon) identifier for each organism recognized by
#' the Biomart service for Ensembl. Internal use.
#'
#' @param org the organism for which to return the identifier.
#' @return A string with the dataset identifier.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' dm3.id <- getDataset("dm3")
#'}
getDataset <- function(org) {
    switch(org,
        hg18 = { return("hsapiens_gene_ensembl") },
        hg19 = { return("hsapiens_gene_ensembl") },
        hg38 = { return("hsapiens_gene_ensembl") },
        mm9 = { return("mmusculus_gene_ensembl") },
        mm10 = { return("mmusculus_gene_ensembl") },
        rn5 = { return("rnorvegicus_gene_ensembl") },
        dm3 = { return("dmelanogaster_gene_ensembl") },
        danrer7 = { return("drerio_gene_ensembl") },
        pantro4 = { return("ptroglodytes_gene_ensembl") },
        susscr3 = { return("sscrofa_gene_ensembl") }
    )
}

#' Annotation downloader helper
#'
#' Returns a vector of chromosomes to maintain after annotation download. Internal
#' use.
#'
#' @param org the organism for which to return the chromosomes. 
#' @return A character vector of chromosomes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg18.chr <- getValidChrs("hg18")
#'}
getValidChrs <- function(org)
{
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        }
    )
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt
#' package in order to fetch the gene annotation for each organism. It has no
#' parameters. Internal use.
#'
#' @param org one of the supported organisms.
#' @return A character vector of Ensembl gene attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' gene.attr <- getGeneAttributes()
#'}
getGeneAttributes <- function(org) {
    if (org %in% c("hg18","mm9"))
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_name",
            "gene_biotype"
        ))
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt
#' package in order to fetch the exon annotation for each organism. It has no
#' parameters. Internal use.
#'
#' @param org one of the supported organisms.
#' @return A character vector of Ensembl exon attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' exon.attr <- getExonAttributes()
#'}
getExonAttributes <- function(org) {
    if (org %in% c("hg18","mm9"))
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else if (org == "hg19")
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

#' Download annotation from UCSC servers, according to organism and source
#'
#' Directly downloads UCSC and RefSeq annotation files from UCSC servers to be
#' used with metaseqR. This functionality is used when the package RMySQL is not
#' available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A data frame with annotation elements.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.file <- getUcscDbl("hg18","gene","ucsc")
#'}
getUcscDbl <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    
    if (!require(RSQLite))
        stop("R package RSQLite is required to use annotation from UCSC!")

    http.base <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/",
        get.ucsc.organism(org),"/database/",sep="")
    table.defs <- getUcscTabledef(org,type,refdb,"fields")
    file.list <- vector("list",length(table.defs))
    names(file.list) <- names(table.defs)
    for (n in names(file.list))
        file.list[[n]] <- paste(http.base,n,".txt.gz",sep="")
        
    # Fill the fields for each table
    drv <- dbDriver("SQLite")
    db.tmp <- tempfile()
    con <- dbConnect(drv,dbname=db.tmp)
    disp("  Retrieving tables for temporary SQLite ",refdb," ",org," ",type,
        " subset database")
    for (n in names(file.list)) {
        disp("    Retrieving table ",n)
        download.file(file.list[[n]],file.path(tempdir(),
            paste(n,".txt.gz",sep="")),quiet=TRUE)
        if (.Platform$OS.type == "unix")
            system(paste("gzip -df",file.path(tempdir(),
                paste(n,".txt.gz",sep=""))))
        else
            unzip(file.path(tempdir(),paste(n,".txt.gz",sep="")))
        sql.df <- read.delim(file.path(tempdir(),paste(n,".txt",sep="")),
            row.names=NULL,header=FALSE,strip.white=TRUE)
        names(sql.df) <- table.defs[[n]]
        dbWriteTable(con,n,sql.df,row.names=FALSE)
    }
    dbDisconnect(con)
    return(db.tmp)
}

#' Get SQLite UCSC table defintions, according to organism and source
#'
#' Creates a list of UCSC Genome Browser database tables and their SQLite
#' definitions with the purpose of creating a temporary SQLite database to be 
#' used used with metaseqR. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A list with SQLite table definitions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.tabledefs <- get.ucsc.tables("hg18","gene","ucsc")
#'}
getUcscTabledef <- function(org,type,refdb="ucsc",what="queries") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    what <- tolower(what[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    checkTextArgs("what",what,c("queries","fields"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    getUcscTblTpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    getUcscTblTpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    getUcscTblTpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    getUcscTblTpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    getUcscTblTpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    getUcscTblTpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                }
            )
        }
    )
}

#' Create SQLite UCSC table template defintions
#'
#' Returns an SQLIte table template defintion, according to  UCSC Genome Browser 
#' database table schemas. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines. Internal use only.
#'
#' @param tab name of UCSC database table.
#' @param what \code{"queries"} for SQLite table definitions or \code{"fields"}
#' for table column names.
#' @return An SQLite table definition.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.table.tmpl <- getUcscTblTpl("knownCanonical")
#'}
getUcscTblTpl <- function(tab,what="queries") {
    if (what=="queries") {
        switch(tab,
            knownCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`knownCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownGene = {
                return(paste(
                    "CREATE TABLE",
                    "`knownGene` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`proteinID` TEXT NOT NULL DEFAULT '',",
                    "`alignID` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            refFlat = {
                return(paste("CREATE TABLE",
                    "`refFlat` (",
                    "`geneName` TEXT NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            knownToEnsembl = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToEnsembl` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            ensemblSource = {
                return(paste(
                    "CREATE TABLE",
                    "`ensemblSource` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`source` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            mgcGenes = {
                return(paste(
                    "CREATE TABLE `mgcGenes` (",
                    "`bin` UNSIGNED INTEGER NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`score` INTEGER DEFAULT NULL,",
                    "`name2` TEXT NOT NULL,",
                    "`cdsStartStat` TEXT NOT NULL,",
                    "`cdsEndStat` TEXT NOT NULL,",
                    "`exonFrames` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            ensemblToGeneName = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToGeneName` (",
                    "`name` TEXT NOT NULL,",
                    "`value` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER unsigned NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            flyBaseGene = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseGene` (",
                    "`bin` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            }
        )
    }
    else if (what=="fields") {
        switch(tab,
            knownCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                "transcript","protein"))
            },
            knownGene = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","proteinID",
                    "alignID"))
            },
            knownToRefSeq = {
                return(c("name","value"))
            },
            refFlat = {
                return(c("geneName","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            knownToEnsembl = {
                return(c("name","value"))
            },
            ensemblSource = {
                return(c("name","source"))
            },
            mgcGenes = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","score",
                    "name2","cdsStartStat","cdsEndStat","exonFrames"
                ))
            },
            ensemblToGeneName = {
                return(c("name","value"))
            },
            flyBaseCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                    "transcript","protein"))
            },
            flyBaseGene = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            flyBaseToRefSeq = {
                return(c("name","value"))
            }
        )
    }
}

#' Return queries for the UCSC Genome Browser database, according to organism and
#' source
#'
#' Returns an SQL query to be used with a connection to the UCSC Genome Browser
#' database and fetch metaseqR supported organism annotations. This query is
#' constructed based on the data source and data type to be returned.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A valid SQL query.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.query <- getUcscQuery("hg18","gene","ucsc")
#'}
getUcscQuery <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is like hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "flyBaseGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`flyBaseCanonical` INNER JOIN `flyBaseGene` ",
                                "ON flyBaseCanonical.transcript=",
                                "flyBaseGene.name INNER JOIN ",
                                "`flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,",
                                "`source` AS `biotype` FROM `refFlat` INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,flyBaseGene.exonStarts AS ",
                                "`start`,flyBaseGene.exonEnds AS `end`,",
                                "`transcript` AS `exon_id`,flyBaseGene.strand ",
                                "AS `strand`,`transcript` AS `gene_id`,",
                                "`geneName` AS `gene_name`,`source` AS ",
                                "`biotype` FROM `flyBaseCanonical` INNER JOIN ",
                                "`flyBaseGene` ON ",
                                "flyBaseCanonical.transcript=flyBaseGene.name ",
                                "INNER JOIN `flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name ",
                                "INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        }
                    )
                }
            )
        }
    )
}

#' Return host, username and password for UCSC Genome Browser database
#'
#' Returns a character vector with a hostname, username and password to connect
#' to the UCSC Genome Browser database to retrieve annotation. Internal use.
#'
#' @return A named character vector.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.creds <- getUcscCredentials()
#'}
getUcscCredentials <- function() {
    return(c(
        host="genome-mysql.cse.ucsc.edu",
        user="genome",
        password=""
    ))
}
