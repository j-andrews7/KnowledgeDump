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

#### Summarization 
The above is nice for generate results for all the significant hits, but it's not a great summary of them. The below will take those results and plot the top X number of significant genesets, ranked by adjusted p-value, for each collection as a plot.

```r
summarize_GSEA <- function(gsea.list, outdir, padj.th = 0.05, top = 75) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  for (i in seq_along(gsea.list)) {
    ct <- names(gsea.list)[i]
    df <- gsea.list[[i]]
    
    df.sub <- df[df$padj < padj.th,]
    
    if (nrow(df.sub) > top) {
      df.sub <- df.sub %>% as_tibble() %>% arrange(padj)
      df.sub <- df.sub[1:top, ]
    }
    
    if (nrow(df.sub) > 0) {
      df.sub <- df.sub %>%
        as_tibble() %>%
        arrange(desc(NES))
      
      p <- ggplot(df.sub, aes(reorder(pathway, -NES), NES)) +
        geom_col(aes(fill=-log10(padj))) + coord_flip() +
        labs(x=NULL, y="Normalized Enrichment Score", 
             title=paste0(ct, " - Top ", top, "\np.adj < ", padj.th)) + 
        theme_bw() + scale_fill_viridis() + ylim(-4,4) + 
        theme(axis.text.y = element_text(size = 6), plot.title = element_text(size = 10)) +
        scale_x_discrete(label = function(x) str_trunc(x, 55))
      
      h <- 2 + (0.07 * nrow(df.sub))
      
      pdf(paste0(outdir, "/", ct, ".padj.", padj.th, ".topbypadj", top, ".revrank.pdf"), width = 7, height = h)
      print(p)
      dev.off()
    }
  }
}

summarize_GSEA(xl.lists, outdir = "./GSEA/RA.v.vehicle")
```

#### Plot Leading Edge Genes
GSEA returns the leading edge genes that are driving the score for a given signature. It can be useful to have a closer look at these genes in the form of boxplots and/or heatmaps.

```r
plot_le <- function(sce, gsea.lists, annot.by, group.by, outdir, use.assay, cells.use, sig.thresh = 0.05, 
					group.by2 = NULL, split.by = NULL, swap.rownames = NULL) {
  
  for (i in seq_along(gsea.lists)) {
    ct <- names(gsea.lists)[i]
    df <- gsea.lists[[i]]
    
    sig.paths <- df$pathway[df$padj < sig.thresh]
    
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    for (p in seq_along(sig.paths)) {
      path.name <- sig.paths[p]
      le <- unlist(df$leadingEdge[df$pathway == path.name])

	  if (nchar(path.name) > 50) {
        path.name <- substr(path.name, 1, 50)
      }

      if (length(le) > 1) {
      
        pdf(paste0(outdir, "/", ct, ".", path.name, ".boxplot.pdf"), width = 5, height = 4)
        
        for (i in use.assay) {
          pl <- dittoPlotVarsAcrossGroups(sce[,cells.use], le, group.by = group.by, 
                                          plots = c("vlnplot", "jitter", "boxplot"), assay = i, sub = i,
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
          print(pl)
          pl <- dittoPlotVarsAcrossGroups(sce[,cells.use], le, group.by = group.by, 
                                          plots = c("vlnplot", "jitter", "boxplot"), 
                                          adjustment = "relative.to.max", assay = i, sub = i,
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
          print(pl)
          pl <- dittoPlotVarsAcrossGroups(sce[,cells.use], le, group.by = group.by, 
                                          plots = c("vlnplot", "jitter", "boxplot"), 
                                          adjustment = "none", assay = i, sub = i,
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
          print(pl)
          
          if (!is.null(group.by2) & !is.null(split.by)) {
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
            print(pl)
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), 
                                            adjustment = "relative.to.max", assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
            print(pl)
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), 
                                            adjustment = "none", assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, swap.rownames = swap.rownames)
            print(pl)
          }
        }
        
        dev.off()
        
        pdf(paste0(outdir, "/", ct, ".", path.name, ".heatmap.pdf"), width = 5, height = 7)
        for (i in use.assay) {
          pl <- dittoHeatmap(sce, le, annot.by = annot.by, cells.use = cells.use, show_colnames = FALSE,
                             breaks = seq(-3, 3, length.out = 51), cluster_rows = FALSE,
                             fontsize_row = 6, cluster_cols = FALSE, assay = i, sub = i, swap.rownames = swap.rownames)
          grid.draw(pl)
          
          pl <- dittoHeatmap(sce, le, annot.by = annot.by, show_colnames = FALSE,
                             breaks = seq(-3, 3, length.out = 51), cluster_rows = FALSE,
                             fontsize_row = 6, cluster_cols = FALSE, assay = i, sub = i, swap.rownames = swap.rownames)
          grid.draw(pl)
        }
        dev.off()
      }
    }
  }
}

plot_le(bulk, xl.lists, annot.by = "Line_Treatment", group.by = "Line_Treatment", group.by2 = "Treatment", 
		split.by = "Line", outdir = "./GSEA/LTC115.v.LTC97_DMSO/leading_edge", use.assay = c("lognorm", "vsd"), 
		cells.use = bulk$Line_Treatment %in% c("LTC115_DMSO", "LTC97_DMSO", "LTC115_Primary", "LTC97_Primary"))
```

### Enrichment/Over-representation Analyses
These are kind of a pain to run for multiple comparisons, etc, so these functions try to ease that pain. Like a glass of water and handful of advil after a rough night.

#### KEGG Enrichment

```r
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("pathview")
library("ReactomePA")

#' @param res.list Named list of DESeq2 results data.frames.
#' @param padj.th Numeric scalar to use as significance threshold for DE genes.
#' @param lfc.th Numeric scalar to use as log fold change threshold for DE genes.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.col Character scalar indicating name of gene ID column for each data.frame in \code{res.list}
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param ... Passed to \code{compareCluster}.
run_enrichKEGG <- function(res.list, padj.th = 0.05, lfc.th = 0, outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.col = "ENSEMBL", id.type = "ENSEMBL", ...) {
  # Do GO enrichment on up/downregulated genes.
  for (r in names(res.list)) {
    df <- res[[r]]
    out <- file.path(outdir, r)
    dir.create(paste0(out, "/KEGG_pathview"), showWarnings = FALSE, recursive = TRUE)
    
    # Strip gene version info if ensembl.
    if (id.type == "ENSEMBL") {
      xx <- strsplit(df[[id.col]], "\\.")
      df[[id.col]] <- unlist(lapply(xx, FUN = function(x) x[1]))
    }
    
    # Get FC values for pathway plots.
    df <- df[!is.na(df$padj),]
    geneList <- bitr(df[[id.col]], fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)
    geneList$FC <- df$log2FoldChange[match(geneList$ENSEMBL, df$ENSEMBL)]
    gl <- geneList$FC
    names(gl) <- geneList$ENTREZID
    gl = sort(gl, decreasing = TRUE)
    
    # Get gene categories.
    genes <- list(upregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange > lfc.th],
                  downregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange < lfc.th],
                  all_de = df[[id.col]][df$padj < padj.th])
    
    genes$upregulated <- bitr(genes$upregulated, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    genes$downregulated <- bitr(genes$downregulated, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    genes$all_de <- bitr(genes$all_de, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    # Remove lowly expressed genes.
    bg <- df[[id.col]][!is.na(df$padj)]
    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    
    ck <- compareCluster(geneCluster = genes, fun = enrichKEGG, keyType = "kegg", universe = bg, ...)
    if (!is.null(ck)) {
      ck <- setReadable(ck, OrgDb = OrgDb, keyType="ENTREZID")
      
      # Term similarities via Jaccard Similarity index.
      ego <- pairwise_termsim(ck)
      
      height = 3 + (0.15 * length(ego@compareClusterResult$Cluster))
      
      pdf(paste0(out, "/KEGG_Enrichments.Top20_perGroup.pdf"), width = 6, height = height)
      p <- dotplot(ego, showCategory = 20, font.size = 9)
      print(p)
      p <- dotplot(ego, size = "count", showCategory = 20, font.size = 9)
      print(p)
      dev.off()
      
      pdf(paste0(out, "/KEGG_Enrichments.termsim.Top20_perGroup.pdf"), width = 9, height = 9)
      p <- emapplot(ego, pie="count", cex_category=0.9, cex_label_category = 0.9, 
                    layout="kk", repel = TRUE, showCategory = 20)
      print(p)
      dev.off()
      
      for (x in ego@compareClusterResult$ID) {
        pname <- ego@compareClusterResult$Description[ego@compareClusterResult$ID == x]
        xx <- tryCatch(
          pathview(gene.data  = gl,
                   pathway.id = x,
                   kegg.dir = paste0(out, "/KEGG_pathview/"),
                   species    = "hsa",
                   out.suffix = pname,
                   limit      = list(gene=3, cpd=1),
                   kegg.native = TRUE),
          error = function(e) {"bah"}
        )
      }

      # This idiotic workaround copies the PNGs to our wanted directory as pathview generates
      # all output in the working directory with no way to alter output location. Pretty dumb.
      old_dir <- "./"
      new_dir <- paste0(out, "/KEGG_pathview/")
      old_files <- list.files(path = old_dir, pattern = "png", full.names = T)
      new_files <- sub(old_dir, new_dir, old_files)
      file.rename(from = old_files, to = new_files)
      
      saveRDS(ego, file = paste0(out, "/enrichKEGG.results.RDS"))
      ego <- as.data.frame(ego)
      write.table(ego, file = paste0(out, "/enrichKEGG.results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

run_enrichKEGG(res)
```

#### Reactome

```r
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("pathview")
library("ReactomePA")

#' @param res.list Named list of DESeq2 results data.frames.
#' @param padj.th Numeric scalar to use as significance threshold for DE genes.
#' @param lfc.th Numeric scalar to use as log fold change threshold for DE genes.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.col Character scalar indicating name of gene ID column for each data.frame in \code{res.list}
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param ... Passed to \code{compareCluster}.
run_enrichPathway <- function(res.list, padj.th = 0.05, lfc.th = 0, outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.col = "ENSEMBL", id.type = "ENSEMBL", ...) {
  # Do GO enrichment on up/downregulated genes.
  for (r in names(res.list)) {
    df <- res[[r]]
    out <- file.path(outdir, r)
    dir.create(paste0(out, "/Reactome_pathways"), showWarnings = FALSE, recursive = TRUE)
    
    # Strip gene version info if ensembl.
    if (id_type == "ENSEMBL") {
      xx <- strsplit(df[[id.col]], "\\.")
      df[[id.col]] <- unlist(lapply(xx, FUN = function(x) x[1]))
    }
    
    # Get FC values for pathway plots.
    df <- df[!is.na(df$padj),]
    geneList <- bitr(df[[id.col]], fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)
    geneList$FC <- df$log2FoldChange[match(geneList$ENSEMBL, df$ENSEMBL)]
    gl <- geneList$FC
    names(gl) <- geneList$ENTREZID
    gl = sort(gl, decreasing = TRUE)
    
    genes <- list(upregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange > lfc.th],
                  downregulated = df[[id.col]][df$padj < padj.th  & df$log2FoldChange < lfc.th],
                  all_de = df[[id.col]][df$padj < padj.th])
    
    genes$upregulated <- bitr(genes$upregulated, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    genes$downregulated <- bitr(genes$downregulated, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    genes$all_de <- bitr(genes$all_de, fromType = id.type, toType = "ENTREZID", 
                      OrgDb = OrgDb)$ENTREZID
    
    # Remove lowly expressed genes.
    bg <- df[[id.col]][!is.na(df$padj)]
    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    
    ck <- compareCluster(geneCluster = genes, fun = enrichPathway, universe = bg, readable = TRUE, ...)
    if (!is.null(ck)) {
      # Term similarities via Jaccard Similarity index.
      ego <- pairwise_termsim(ck)
      
      # Adjust plot height based on number of terms.
      height = 4 + (0.1 * length(ego@compareClusterResult$Cluster))
      
      pdf(paste0(out, "/Reactome_Enrichments.Top20_perGroup.pdf"), width = 6, height = height)
      p <- dotplot(ego, showCategory = 20, font.size = 9)
      print(p)
      p <- dotplot(ego, size = "count", showCategory = 20, font.size = 9)
      print(p)
      dev.off()
      
      pdf(paste0(out, "/Reactome_Enrichments.termsim.Top20_perGroup.pdf"), width = 9, height = 9)
      p <- emapplot(ego, pie="count", cex_category=0.9, cex_label_category = 0.9, 
                    layout="kk", repel = TRUE, showCategory = 20)
      print(p)
      dev.off()
      
      # Remove duplicate IDs.
      gl <- gl[unique(names(gl))]
      
      # Network plot for each pathway with FC values.
      for (x in ego@compareClusterResult$Description) {
        # Some pathways have backslashes, which will break file creation.
        x_out <- str_replace_all(x, "/", "_")
        
        # For when it inevitably wants to crash due to not finding a pathway name or such.
        tryCatch(
          {
            pdf(paste0(out, "/Reactome_pathways/", x_out, ".pdf"), width = 11, height = 11)
            p <- viewPathway(x, readable = TRUE, foldChange = gl)
            vals <- p$data$color[!is.na(p$data$color)]
            l <- max(abs(as.numeric(vals)))
            p <- p + scale_color_gradient2(limits = c(-l,l), mid = "grey90", 
                                           high = "red", low = "navyblue")
            print(p)
            dev.off()
            dev.off()
            dev.off()
          },
          error = function(e) {"bah"}
        )
      }
      
      saveRDS(ego, file = paste0(out, "/enrichPathway.reactome.RDS"))
      ego <- as.data.frame(ego)
      write.table(ego, file = paste0(out, "/enrichPathway.reactome.results.txt"), 
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}
```

### CNV Calling from Methylation Array
This spits out typical genome-wide CNV plots, segmentation files, bins, and IGV tracks from Illumina methylation arrays. Users can add details regions for labels if they'd like. When mixing both 450k and EPIC arrays, set `array_type = "overlap"`.
```r
library("minfi")
library("conumee")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

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
    anno@probes <- anno@probes[names(anno@probes) %in% 
      names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]
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

