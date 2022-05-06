---
aliases: [code snippets and functions]
---

# Useful Code Snippets and Functions
This serves as a collection of code snippets and functions for common data science, bioinformatics, data cleaning/munging, and data visualization tasks. I got sick of hunting through old notebooks and scripts to find that one thing I wrote a year ago.

# R

## Gene/Peak Annotations & Conversions

### Convert Ensembl Gene IDs to Symbols (or vice versa)
```r
# Mouse
library(ensembldb)
library(EnsDb.Mmusculus.v79)

ens.ids <- c("ENSMUSG00000025907")

symbs <- mapIds(EnsDb.Mmusculus.v79, keys = ens.ids, keytype = "GENEID", column = "SYMBOL")
symbs
```

## Viz

### 3D tSNE/UMAP/PCA

```r
library(plotly)
library(SingleCellExperiment)
library(dittoSeq)

#' @param sce SingleCellExperiment object.
#' @param dimred Character scalar indicating the name of the dimensionality reduction to plot.
#' @param color.by Character scalar indicating colData column to color points by.
#' @param shape.by Character scalar indicating colData column to shape points by.
#' @param hover.info Character scalar or vector indicating colData column(s) to display when points are hovered.
#' @param pt.size Numeric scalar indicating point size.
plot3Ddim <- function(sce, dimred, color.by = NULL, shape.by = NULL, hover.info = NULL, pt.size = 3) {
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
      hov[[i]] <- paste0("</br><b>",hover.info[i], ":</b> ", colData(sce)[[hover.info[i]]])
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

fig <- plot3Ddim(new.sce, "UMAP_m.dist0.3_n.neigh10", color.by = "Group", hover.info = c("CellType", "Group"))

# Self as self-contained html if wanted.
saveWidget(jqui_resizable(fig), "./QC/UMAP.3D.m.dist0.3_n.neigh10.Group.html")
```

## Single Cell RNA-seq
### dimReduc Sweep
To get lots of dimensionality reductions with differing parameters.

```r
library(SingleCellExperiment)
library(scater)

umap_sweep <- function(sce, dim_reduc, min_dist = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3), n_neighbors = c(10, 15, 20, 30, 40, 50)) {
  for (d in min_dist) {
    for (n in n_neighbors) {
      sce <- runUMAP(sce, n_neighbors = n, min_dist = d, name = paste0("UMAP_m.dist", d, "_n.neigh", n), dimred = dim_reduc, ncomponents = 2, BPPARAM = SnowParam(6))
    }
  }
  
  return(sce)
}

sce <- umap_sweep(sce, dim_reduc = "corrected")
```

### cluster Sweep
To get lots of different clustering sets with different methods/parameters.

```r
library(SingleCellExperiment)
library(scater)
library(bluster)
library(dittoSeq)

out <- clusterSweep(reducedDim(sce, "corrected"), 
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

# Python


# Bash

### Extract Any Type of Compressed File
I just add this to my `.bashrc` file.

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
