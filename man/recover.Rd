\name{recover}
\alias{recover}
\title{Create genomic signal profiles in predefined or 
custom areas using short sequence reads}
\usage{
    recover <- function(
        input,
        design = NULL,
        region = c("genebody", "tss", "tes", "custom"),
        type = c("chipseq", "rnaseq"),
        genome = c("hg18", "hg19", "hg38", "mm9" ,"mm10",
            "rn5", "dm3", "danrer7", "pantro4", "susscr3"),
        refdb = c("ensembl", "ucsc", "refseq"),
        flank = c(2000, 2000),
        orderBy = list(
            what = c("none", "suma", "sumn", "maxa",
                "maxn", "hcn"),
            order = c("descending", "ascending")
        ),
        binParams = list(
            flankBinSize = 0,
            regionBinSize = 0,
            sumStat = c("mean", "median"),
            smooth = TRUE,
            forceHeatmapBinning = TRUE,
            forcedBinSize = c(50, 200)
        ),
        selector = list(
            id = NULL,
            biotype = getBiotypes(genome),
            exonType = NULL
        ),
        preprocessParams = list(
            normalize = c("none", "linear", "downsample",
                "sampleto"),
            sampleTo = 1e+6,
            spliceAction = c("keep","remove","split"),
            spliceRemoveQ = 0.75,
            seed = 42
        ),
        plotParams = list(
            profile = TRUE,
            heatmap = TRUE,
            signalScale = c("natural", "log2"),
            heatmapScale = c("each", "common"),
            device = c("x11", "png", "jpg", "tiff", "bmp",
                "pdf", "ps"),
            outputDir = ".",
            outputBase = NULL
        ),
        saveParams = list(
            ranges = TRUE,
            coverage = TRUE,
            profile = FALSE
        ),
        complexHeatmapParams = list(
            main = list(
                cluster_rows = ifelse(orderBy$what == "hc",
                    TRUE, FALSE),
                cluster_columns = FALSE,
                column_title_gp = gpar(fontsize = 12, 
                    font = 2),
                show_row_names = FALSE,
                show_column_names = FALSE,
                heatmap_legend_param = list(
                    color_bar = "continuous"
                )
            ),
            group = list(
                cluster_rows = ifelse(orderBy$what == "hc",
                    TRUE, FALSE),
                cluster_columns = FALSE,
                column_title_gp = gpar(fontsize = 12,
                    font = 2),
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 10, 
                    font = 2),
                gap = unit(5,"mm"),
                heatmap_legend_param = list(
                    color_bar = "continuous"
                )
            )
        ),
        kmParams = list(
            k = 0,
            nstart = 20,
            algorithm = c("Hartigan-Wong", "Lloyd",
                "Forgy", "MacQueen"),
            iterMax = 20,
            reference = NULL,
            seed = 42
        ),
        strandedParams = list(
            strand = NULL,
            ignoreStrand = TRUE
        ),
        bamParams = NULL,
        onTheFly = FALSE,
        localDbHome = file.path(path.expand("~"), ".recover"),
        rc = NULL
    )
}
\arguments{
    \item{input}{}

    \item{design}{}

    \item{region}{}
    
    \item{type}{}
    
    \item{genome}{}
    
    \item{refdb}{}

    \item{flank}{}

    \item{orderBy}{}
    
    \item{binParams}{}
    
    \item{selector}{}
    
    \item{preprocessParams}{}

    \item{plotParams}{}

    \item{saveParams}{}
    
    \item{complexHeatmapParams}{}
    
    \item{kmParams}{}
    
    \item{strandedParams}{}

    \item{bamParams}{BAM file read parameters passed to 
    \code{\link{BamFile}. See the related function.
    Currently this is not used.}

    \item{onTheFly}{Read short reads directly from BAM files
    when input contains paths to BAM files. In this case the
    storage of short reads in the output list object as a 
    GRanges object is not possible.}
    
    \item{localDbHome}{}
    
    \item{rc}{}
}
\value{
    a named list with 
}
\description{
    This function performs calculates and plots signal
    profiles created from short sequence reads derived
    from Next Generation Sequencing technologies. The
    profiles provided are either sumarized curve profiles
    or heatmap profiles. Currently, \code{recover} 
    supports genomic profile plots for reads derived from
    ChIP-Seq and RNA-Seq experiments.
}
\examples{
\donttest{
}
\author{
    Panagiotis Moulos
}

