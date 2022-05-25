---
aliases: [code snippets and functions]
---

# Useful Code Snippets and Functions
This serves as a collection of code snippets and functions for common data science, bioinformatics, data cleaning/munging, and data visualization tasks. I got sick of hunting through old notebooks and scripts to find that one thing I wrote a year ago.

## R

### Gene/Peak Annotations & Conversions

#### Convert Ensembl Gene IDs to Symbols (or vice versa)
There's like 45 ways to do this, but these are pretty easy with a decent recovery rate.

=== "Using mygene"

    ```r
    library(mygene)
	
	# Species can be set with `query` and `queryMany`. Symbol to ensembl.
	# Set `fields = "all"` to get all the info available.
	df <- queryMany(c("Cdk2", "Cdk3"), fields = c("ensembl.gene", "entrezgene"), species = "mouse", size = 1)
    df$ensembl.gene
    ```

=== "Using ensembldb"
	```r
	library(ensembldb)
	library(EnsDb.Mmusculus.v79)
	
	ens.ids <- c("ENSMUSG00000025907")
	
	symbs <- mapIds(EnsDb.Mmusculus.v79, keys = ens.ids, keytype = "GENEID", column = "SYMBOL")
	symbs
	```

### Gene Summaries
Refseq description summaries. Can use almost any ID for these, does a decent job figuring it out.

```r
library(mygene)

# One gene.
gene <- query("CDK2", fields = "all", size = 1)
gene$hits$summary


# Multiple genes. Returns a dataframe. Set `fields = "all"` to get all the info available.
df <- queryMany(c(1017, 1018, "ENSG00000148795", "LAIR1"), 
               fields = c("symbol", "name", "taxid", "entrezgene", "summary"), size = 1)
df$summary
```

### Viz

#### 3D tSNE/UMAP/PCA
Occasionally, 3D dimensionality reduction plots can be kind of useful for exploratory purposes.

```r
library(plotly)
library(SingleCellExperiment)
library(dittoSeq)

#' @param sce SingleCellExperiment object.
#' @param dimred Character scalar indicating the name of the dimensionality reduction to plot.
#' @param color.by Character scalar indicating colData column to color points by.
#' @param shape.by Character scalar indicating colData column to shape points by.
#' @param hover.info Character scalar or vector indicating colData column(s) 
#'   to display when points are hovered.
#' @param pt.size Numeric scalar indicating point size.
plot3Ddim <- function(sce, dimred, color.by = NULL, shape.by = NULL, 
					  hover.info = NULL, pt.size = 3) {
					  
  dimmy <- as.data.frame(reducedDim(sce, dimred))
  names(dimmy) <- paste0(dimred, "_", 1:3)
  
  pl.shapes <- NULL
  pl.cols <- NULL
  pl.col <- "black"
  
  if (!is.null(color.by)) {
    pl.cols <- colData(sce)[,color.by, drop = TRUE]
    
    if (is.factor(pl.cols)) {
      pl.cols <- droplevels(pl.cols)
    }
  }
  
  if (!is.null(shape.by)) {
    pl.shapes <- colData(sce)[,shape.by, drop = TRUE]
    
    if (is.factor(pl.cols)) {
      pl.shapes <- droplevels(pl.shapes)
    }
  }
  
  pl.col <- dittoColors()[seq_along(unique(colData(sce)[,color.by, drop = TRUE]))]
  
  hov.text <- NULL
  
  if (!is.null(hover.info)) {
    hov <- list()
    for (i in seq_along(hover.info)) {
      hov[[i]] <- paste0("</br><b>",hover.info[i], ":</b> ", 
                         colData(sce)[[hover.info[i]]])
    }
  
    hov.text <- do.call(paste0, hov)
  }
  
  # Generate plot.
  fig <- plot_ly(dimmy, x = as.formula(paste0("~", dimred, "_1")),
          y = as.formula(paste0("~", dimred, "_2")),
          z = as.formula(paste0("~", dimred, "_3")),
          type = "scatter3d",
          mode = "markers",
          color = pl.cols,
          colors = pl.col,
          size = pt.size,
          symbol = pl.shapes,
          symbols = c("circle", "square", "diamond", "cross", "diamond-open",
                      "circle-open", "square-open", "x"),
          text = hov.text,
          hoverinfo = "text") %>%
    layout(scene = list(
      xaxis = list(title = paste0(dimred, "_1")),
      yaxis = list(title = paste0(dimred, "_2")),
      zaxis = list(title = paste0(dimred, "_3")),
      camera = list(eye = list(x=1.5, y = 1.8, z = 0.4)))) %>%
    toWebGL()
  
  return(fig)
}

fig <- plot3Ddim(new.sce, "UMAP_m.dist0.3_n.neigh10", 
				 color.by = "Group", hover.info = c("CellType", "Group"))

# Self as self-contained html if wanted.
saveWidget(jqui_resizable(fig), "./QC/UMAP.3D.m.dist0.3_n.neigh10.Group.html")
```

### Single Cell RNA-seq
#### dimReduc Sweep
To get lots of dimensionality reductions with differing parameters.

```r
library(SingleCellExperiment)
library(scater)
library(BiocParallel)

#' @param sce SingleCellExperiment object.
#' @param dimred Character scalar indicating the name of the dimensionality reduction use as input.
#' @param min_dist Numeric vector indicating parameters to sweep for min_dist UMAP parameter.
#' @param n_neighbors Numeric vector indicating parameters to sweep for n_neighbors UMAP parameter.
umap_sweep <- function(sce, dim_reduc, 
					   min_dist = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3), 
					   n_neighbors = c(10, 15, 20, 30, 40, 50)) {
					   
  for (d in min_dist) {
    for (n in n_neighbors) {
      sce <- runUMAP(sce, n_neighbors = n, min_dist = d, 
                     name = paste0("UMAP_m.dist", d, "_n.neigh", n), 
	                 dimred = dim_reduc, ncomponents = 2, BPPARAM = SnowParam(6))
    }
  }
  
  return(sce)
}

sce <- umap_sweep(sce, dim_reduc = "PCA")
```

#### cluster Sweep
To get lots of different clustering sets with different methods/parameters.

```r
library(SingleCellExperiment)
library(scater)
library(bluster)
library(dittoSeq)

out <- clusterSweep(reducedDim(sce, "PCA"), 
    NNGraphParam(), 
    k=as.integer(c(10, 15, 20, 25, 30, 35, 40)),
    cluster.fun=c("louvain", "walktrap"))

# Cluster metrics
df <- as.data.frame(out$parameters)
df$num.clusters <- vapply(as.list(out$clusters), function(cluster) { 
    length(unique(cluster))
}, 0L)

all.sil <- lapply(as.list(out$clusters), function(cluster) {
    sil <- approxSilhouette(reducedDim(new.sce), cluster)
    mean(sil$width)
})
df$silhouette <- unlist(all.sil)

all.wcss <- lapply(as.list(out$clusters), function(cluster) {
    sum(clusterRMSD(reducedDim(new.sce), cluster, sum=TRUE), na.rm=TRUE)
})
df$wcss <- unlist(all.wcss)

# Plotting cluster metrics
pdf("./QC/clustering.corr.pdf", height = 4, width = 12)
gridExtra::grid.arrange(
    ggplot(df, aes(x=k, y=num.clusters, group=cluster.fun, color=cluster.fun)) + 
        geom_line(lwd=2) + scale_y_log10(),
    ggplot(df, aes(x=k, y=silhouette, group=cluster.fun, color=cluster.fun)) + 
        geom_line(lwd=2),
    ggplot(df, aes(x=k, y=wcss, group=cluster.fun, color=cluster.fun)) + 
        geom_line(lwd=2),
    ncol=3
)
dev.off()

# Add to SCE object.
celldata <- colData(sce)
colData(sce) <- cbind(celldata, out$clusters)
```


### GSEA
#### High-throughput Functions with Plotting
This will run through a bunch of named lists containing genesets and run GSEA on them, yielding both text output for each geneset collection and plots. Two variants, one specific for `msigdbr` genesets (`runGSEA`), one for custom genesets (`runCustomGSEA`).

I'll write a tutorial for doing this end-to-end eventually.

```r
library("fgsea")
library("msigdbr")
library("dplyr")
library("gridExtra")
library("BiocParallel")

#' Run GSEA via fgsea.
#'
#' @param msigs A set of gene sets as returned by \code{msigdbr}.
#' @param ranked.genes A vector of gene identifiers ranked by a test statistic, effect size, etc.
#' @param outdir Character scalar indicating the output directory.
#' @param outprefix Character scalar indicating the output prefix for files.
#' @param xlsx A named list containing results from previous runs. Allows for
#'   more brainless looping.
#' @param cats A character scalar or vector containing MSigDB categories
#'   to use from full set of gene sets.
#' @param subcats A character scalar or vector for subcategories within each category to limit to.
#'   If provided, should be the same length as \code{cats}.
#' @param ... Additional arguments passed to \code{fgsea}.
#' @return xlsx, a named list of GSEA results for each (sub)category.
#'
runGSEA <- function(msigs, ranked.genes, outdir, outprefix, 
                    xlsx = NULL, cats = "H", subcats = NULL, ...) {
  if (length(cats) != length(subcats)) {
    stop("cats and subcats must be of equal length")
  }
  
  collapsedPathways <- NULL
  
  for (i in seq_along(cats)) {
    # Parse out the category and subcategory, set output prefixes.
    subcat <- NULL
    categ <- cats[i]
    if (!is.null(subcats) & subcats[i] != "") {
      subcat <- subcats[i]
      sigs <- msigs %>% dplyr::filter(gs_cat == categ & gs_subcat == subcat)
      outpre <- paste0(outdir, "/", outprefix, ".", categ, ".", subcat)
      outpre <- gsub(":", ".", outpre)
      xlname <- paste0(outprefix, ".", categ, ".", subcat)
      xlname <- gsub(":", ".", xlname)
    } else {
      sigs <- msigs %>% dplyr::filter(gs_cat == categ)
      outpre <- paste0(outdir, "/", outprefix, ".", categ)
      xlname <- paste0(outprefix, ".", categ)
    }
    
    # Convert the msigdbr dataframe to named lists containing the gene sets in the set category.
    sigs <- sigs %>% split(x = .$gene_symbol, f = .$gs_name)
    fgseaRes <- fgsea(pathways = sigs, 
                  stats    = ranked.genes,
                  eps      = 1e-100,
                  minSize  = 15,
                  maxSize  = 1000,
                  ...)
    fgseaRes <- fgseaRes[order(padj),]
    
    # Save full results.
    fwrite(fgseaRes, file=paste0(outpre, ".fgseaRes.txt"), sep="\t", sep2=c("", " ", ""))
      
    # Figures
    fsig <- fgseaRes$pathway[fgseaRes$padj < 0.05]
    plots <- list()
    for (f in seq_along(fsig)) {
      pathw <- fsig[f]
      if (!is.na(pathw)) {
        # Adjust corner that stats will be plotted in based on swoop shape.
        if (!is.na(fgseaRes$NES[f]) & fgseaRes$NES[f] < 0) {
          xinf <- -Inf
          yinf <- -Inf
        } else {
          xinf <- Inf
          yinf <- Inf
        }
        
        # For those really long titles.
        tt <- pathw
        if (nchar(tt) > 40) {
          stri_sub(tt, 48, 47) <- "\n"
        } 
        
        p <- plotEnrichment(sigs[[pathw]],
                ranked.genes) + labs(title = tt) + 
        theme(plot.title = element_text(size=6))

        # Add stats to plot.
        p <- p + annotate("text", xinf, yinf,
                          label = paste0("p.val = ",
                                         formatC(fgseaRes$pval[fgseaRes$pathway == pathw],
                                                 format = "e", digits = 2),
                                         "\np.adj = ",
                                         formatC(fgseaRes$padj[fgseaRes$pathway == pathw],
                                               format = "e", digits = 2),
                                         "\nNES = ",
                                         round(fgseaRes$NES[fgseaRes$pathway == pathw], digits = 2)),
                        vjust = "inward", hjust = "inward", size = 3)
        plots[[f]] <- p
      }
    }
    
    if (length(plots) > 0) {
      pdf(paste0(outpre, ".Pathways.padj0.05.Swoops.pdf"), height = 10, width = 20)
      # Calculate how many pages to print assuming max 24 plots per page.
      pages <- ceiling(length(plots)/24)
      # Print each page.
      for (i in 1:pages) {
        end <- i * 24
        start <- end - 23
        if (end > length(plots)) {
          end <- length(plots)
        }
        grid.arrange(grobs = plots[start:end], nrow = 4, ncol = 6)
      }
      dev.off()
    }
    
    # Plot top 10 pos/neg enriched pathways in table-ish format plot.
    if (length(fsig) > 0) {
      pdf(paste0(outpre, ".Top10Pathways.padj0.05.pdf"), width = 12)
      topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.05][head(order(pval), n=10), pathway]
      topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.05][head(order(pval), n=10), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      plotGseaTable(sigs[topPathways], ranked.genes, fgseaRes, 
              gseaParam=0.5)
      dev.off()
    }
    
    # Add results to named list.
    if (!is.null(xlsx)) {
      xlsx[[xlname]] <- fgseaRes
    }
  }
  if (!is.null(xlsx)) {
      return(xlsx)
  }
}


runCustomGSEA <- function(sigs, ranked.genes, outdir, 
                          outprefix, xlsx = NULL, ...) {
  # Basically the same function as above except 'sigs' is just
  # a named list of gene sets.
  
  outpre <- paste0(outdir, "/", outprefix, ".c")
  xlname <- paste0(outprefix, ".c")

  fgseaRes <- fgsea(pathways = sigs, 
                stats    = ranked.genes,
                eps      = 1e-100,
                minSize  = 15,
                maxSize  = 3000,
                ...)
  fgseaRes <- fgseaRes[order(padj),]
  
  # Save full results.
  fwrite(fgseaRes, file=paste0(outpre, ".fgseaRes.txt"), 
         sep="\t", sep2=c("", " ", ""))
    
  # Figures
  fsig <- fgseaRes$pathway
  plots <- list()
  for (f in seq_along(fsig)) {
    pathw <- fsig[f]
    if (!is.na(pathw)) {
      
      # reposition annotation depending on curve shape.
      if (!is.na(fgseaRes$NES[f]) & fgseaRes$NES[f] < 0) {
        xinf <- -Inf
        yinf <- -Inf
      } else {
        xinf <- Inf
        yinf <- Inf
      }
      
      # For those really long titles.
      tt <- pathw
      if (nchar(tt) > 40) {
        stri_sub(tt, 48, 47) <- "\n"
      } 
      
      p <- plotEnrichment(sigs[[pathw]],
              ranked.genes) + labs(title=tt) + 
        theme(plot.title = element_text(size=7))
      
      # Add stats to plot.
      p <- p + annotate("text", xinf, yinf, 
                        label = paste0("p.val = ", 
                                       formatC(fgseaRes$pval[f], format = "e", digits = 2), 
                                       "\np.adj = ", 
                                       formatC(fgseaRes$padj[f], format = "e", digits = 2), 
                                       "\nNES = ",
                                       round(fgseaRes$NES[f], digits = 2)),
                        vjust = "inward", hjust = "inward", size = 3)
      
      plots[[f]] <- p
    }
  }
  
  pdf(paste0(outpre, ".Pathways.padj0.05.Swoops.pdf"), height = 10, width = 20)
  # Calculate how many pages to print assuming max 24 plots per page.
  pages <- ceiling(length(plots)/24)
  # Print each page.
  for (i in 1:pages) {
    end <- i * 24
    start <- end - 23
    if (end > length(plots)) {
      end <- length(plots)
    }
    grid.arrange(grobs = plots[start:end], nrow = 4, ncol = 6)
  }
  dev.off()
  
  pdf(paste0(outpre, ".Top10Pathways.padj0.05.pdf"), width = 12)
  topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.05][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.05][head(order(pval), n=10), pathway]
  topPathways <- unique(c(topPathwaysUp, rev(topPathwaysDown)))
  plotGseaTable(sigs[topPathways], ranked.genes, fgseaRes, 
          gseaParam=0.5)
  dev.off()
  
  if (!is.null(xlsx)) {
    xlsx[[xlname]] <- fgseaRes
  }
  
  if (!is.null(xlsx)) {
    return(xlsx)
  }
}
```

### CNV Calling from Methylation Array
```r
library("minfi")
library("conumee")

#' @param meta Character scalar for theath to samplesheet with sample metadata, each row containing a sample.
#' @param basedir Character scalar for the base directory containing IDATs.
#' @param controls Character scalar or vector.
#' @param outdir Character scalar or vector.
#' @param exclude_regions GRanges object containing regions to exclude from the CN plots.
#' @param detail_regions GRanges object containing regions to label.
#' @param array_type Character scalar indicating more array type. Options are "450k", "EPIC", or "overlap" for
#'   datasets with mixed arrays.
#' @param idat_cols Character scalar or vector containing column names for sample identifier columns to paste together.
#'   This should be the IDAT ID.
#' @param name_col Character scalar for column containing sample names.
run_conumee_CNV <- function(meta, basedir, controls, outdir, exclude_regions = NULL, detail_regions = NULL, 
                            array_type = "450k", idat_cols = c("Sentrix_ID", "Sentrix_Position"),
                            name_col = "Sample") {
  ## Load data.
  meta <- read.csv(meta)
  meta$Basename <- file.path(basedir, apply(meta[,idat_cols, drop = FALSE], MARGIN = 1, FUN = paste0, collapse = ""))
  samps <- read.metharray.exp(targets = meta, force = TRUE)
  
  samps <- preprocessNoob(samps)
  
  anno <- CNV.create_anno(array_type = array_type, exclude_regions = exclude_regions, detail_regions = detail_regions)
  
  # This is a bugfix for EPIC arrays and the probes in the annotations by default being screwed up. 
  if (array_type %in% c("EPIC", "overlap")) {
    anno@probes <- anno@probes[names(anno@probes) %in% names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]
  }
  
  ## CNV Calling
  cn.data <- CNV.load(samps)

  sammies <- unlist(pData(samps)[name_col])
  sammies <- sammies[!sammies %in% controls]
  
  ## Plots & Tables
  pdf(paste0(outdir, "/CNVplots.conumee.pdf"), height = 9, width = 18)
  for (s in sammies) {
    s.data <- rownames(pData(samps))[unlist(pData(samps)[name_col]) == s]
    c.data <- rownames(pData(samps))[unlist(pData(samps)[name_col]) %in% controls]
    x <- CNV.fit(cn.data[s.data], cn.data[c.data], anno)
    
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    
    CNV.genomeplot(x, main = s)
    
    CNV.write(x, what = "segments", file = paste0(outdir, "/", s, ".CNVsegments.seg"))  
    CNV.write(x, what = "bins", file = paste0(outdir, "/", s, ".CNVbins.igv"))
    CNV.write(x, what = "detail", file = paste0(outdir, "/", s, ".CNVdetail.txt"))
    CNV.write(x, what = "probes", file = paste0(outdir, "/", s, ".CNVprobes.igv"))
  }
  dev.off()

}

## Detail regions can be made from a BED file if wanted, see the example data for format.
data(exclude_regions)
data(detail_regions)

run_conumee_CNV(meta = "SampleMap.csv", 
                array_type = "EPIC",
                basedir = "./", 
                controls = c("Ctrl1", "Ctrl2", "Ctrl3"), 
                outdir = "./cnv",
                exclude_regions = exclude_regions,
                detail_regions = detail_regions,
                idat_cols = "IDAT",
                name_col = "Sample")
```

## Python
:snake: This is content.

## Bash/Unix Tools
### Useful `.bashrc` Functions, Aliases, etc.
Stuff to cram into your `.bashrc` to make your life generally easier.

#### Extract Any Type of Compressed File

```bash
function extract {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract ."
 else
if [ -f $1 ] ; then
        # NAME=${1%.*}
        # mkdir $NAME && cd $NAME
        case $1 in
          *.tar.bz2) tar xvjf $1 ;;
          *.tar.gz) tar xvzf $1 ;;
          *.tar.xz) tar xvJf $1 ;;
          *.lzma) unlzma $1 ;;
          *.bz2) bunzip2 $1 ;;
          *.rar) unrar x -ad $1 ;;
          *.gz) gunzip $1 ;;
          *.tar) tar xvf $1 ;;
          *.tbz2) tar xvjf $1 ;;
          *.tgz) tar xvzf $1 ;;
          *.zip) unzip $1 ;;
          *.Z) uncompress $1 ;;
          *.7z) 7z x $1 ;;
          *.xz) unxz $1 ;;
          *.exe) cabextract $1 ;;
          *) echo "extract: '$1' - unknown archive method" ;;
        esac
else
echo "$1 - file does not exist"
    fi
fi
}
```

