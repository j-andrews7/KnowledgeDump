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

### Convert human to mouse gene orthologs

=== "Using babelgene"
    ```r
    library(babelgene)
	orthologs(genes, species, human = TRUE, min_support = 3, top = TRUE)
	```

=== "Using biomart"
    ```r
    library(biomart)
    convertHumanGeneList <- function(x) {
	   human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , 
		  mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	  humanx <- unique(genesV2[, 2])
	  return(humanx)
	}
	genes <- convertMouseGeneList(humGenes)
	```

=== "Using JAX homolog list"
	```r
	library(dplyr)
	mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
	
	convert_mouse_to_human <- function(gene_list) {
	
	  output = c()
	
	  for(gene in gene_list){
	    class_key = (mouse_human_genes %>% 
		    filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
	    if(!identical(class_key, integer(0))) {
	      human_genes = (mouse_human_genes %>% 
		      filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
	      for(human_gene in human_genes){
	        output = append(output,human_gene)
	      }
	    }
	  }
	
	  return (output)
	}
	
	convert_human_to_mouse <- function(gene_list){
	
	  output = c()
	
	  for(gene in gene_list){
	    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
	    if(!identical(class_key, integer(0)) ){
	      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
	      for(mouse_gene in mouse_genes){
	        output = append(output, mouse_gene)
	      }
	    }
	  }
	
	  return (output)
	}
	```


Biomart goes down pretty often, so the `babelgene` option is more reliable.

### Viz

#### 3D tSNE/UMAP/PCA
Occasionally, 3D dimensionality reduction plots can be kind of useful for exploratory purposes.

```r
library(plotly)
library(SingleCellExperiment)
library(dittoSeq)
library(htmlwidgets)
library(shinyjqui)


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

#### Gene Expression Over Time (or Between Groups)

Useful for time-series RNA-seq and such.

```r
library(dittoSeq)
library(ggplot2)

# Single gene. Add "method = 'lm'" to 'geom_smooth' for straight-line trend.
dittoPlot(dds, "Clu", group.by = "Broad_Group", color.by = "H3_Status", split.by = "Location", 
		  assay = "lognorm", swap.rownames = "SYMBOL", adjustment = NULL, 
		  plots = c("boxplot", "jitter"), boxplot.width = 0.6, boxplot.lineweight = 0.5, 
		  ylab = "log2(normalized counts + 1)") + 
		  geom_smooth(aes(group = color, color = color), se = FALSE, linewidth = 1.5) + 
		  scale_color_manual(values = Darken(dittoColors()[1:2]))

# Group of genes.
dittoPlotVarsAcrossGroups(dds, c("Cdk2", "Jun", "Fos", "Fcmr"), group.by = "Broad_Group", 
						  color.by = "H3_Status", split.by = "Location", assay = "lognorm", 
						  swap.rownames = "SYMBOL", adjustment = NULL, plots = c("boxplot", "jitter"), 
						  boxplot.width = 0.6, main = "Jessa.Pons.RGC_and_progenitor", 
						  boxplot.lineweight = 0.5) + 
						  geom_smooth(aes(group = color, color = color), se = FALSE, linewidth = 1.5) + 
						  scale_color_manual(values = Darken(dittoColors()[1:2]))
```

#### plotly Subplot Orientation Mirroring

For when you want multiple subplots to mirror each other in terms of orientation. Note that `scene` must be set properly in each `plot_ly` call.

```r
library(plotly)

pal <- c("#18B803", "#138901", "#b3b3b3", "#5A5A5A")
pal <- setNames(pal, c("s1", "s2", "s3", "s4"))

fig1 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~sample, mode = "markers", marker = list(size = 3), scene = "scene1", colors = pal) %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "sample", showarrow = F, 
xref='paper', yref='paper'))

fig2 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~Vim, mode = "markers", marker = list(size = 3), scene = "scene2") %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "Vim", showarrow = F, 
xref='paper', yref='paper'))

fig3 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~Mbp, mode = "markers", marker = list(size = 3), scene = "scene3") %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "Mbp", showarrow = F, 
xref='paper', yref='paper'))

fig4 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~Pdgfra, mode = "markers", marker = list(size = 3), scene = "scene4") %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "Pdgfra", showarrow = F, 
xref='paper', yref='paper'))

fig5 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~Plp1, mode = "markers", marker = list(size = 3), scene = "scene5") %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "Plp1", showarrow = F, 
xref='paper', yref='paper'))

fig6 <- mdf.rpca %>% plot_ly(x = ~DC1, y = ~DC2, z = ~DC3, color = ~Fyn, mode = "markers", marker = list(size = 3), scene = "scene6") %>% 
     layout(annotations = list(x = 0.2 , y = 1.05, text = "Fyn", showarrow = F, 
xref='paper', yref='paper'))

main_plot <- subplot(fig1, fig2, fig3, fig4, fig5, fig6, nrows = 2, margin = 0.06) %>% 
  layout(
         scene  = list(domain = list(x = c(0, 0.33), y = c(0.5, 1)), aspectmode = "cube"), 
         scene2 = list(domain = list(x = c(0.33, 0.66), y = c(0.5, 1)), aspectmode = "cube"),
         scene3 = list(domain = list(x = c(0.66, 1), y = c(0.5, 1)), aspectmode = "cube"),
         scene4 = list(domain = list(x = c(0, 0.33), y = c(0, 0.5)), aspectmode = "cube"),
         scene5 = list(domain = list(x = c(0.33, 0.66), y = c(0, 0.5)), aspectmode = "cube"),
         scene6 = list(domain = list(x = c(0.66, 1), y = c(0, 0.5)), aspectmode = "cube")
  )

main_plot <- main_plot %>% 
  htmlwidgets::onRender(
    "function(x, el) {
      x.on('plotly_relayout', function(d) {
        const camera = Object.keys(d).filter((key) => /\\.camera$/.test(key));
        if (camera.length) {
          const scenes = Object.keys(x.layout).filter((key) => /^scene\\d*/.test(key));
          const new_layout = {};
          scenes.forEach(key => {
            new_layout[key] = {...x.layout[key], camera: {...d[camera]}};
          });
          Plotly.relayout(x, new_layout);
        }
      });
    }")
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

#### Downsample an SCE

For testing stuff on smaller numbers of cells, etc.

```r
#' Downsample a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object.
#' @param ncells Number of cells to downsample to.
#' @return Downsampled SingleCellExperiment object.
downsampleSCE <- function(sce, ncells) {
	keep <- sample(seq(1,ncol(sce),by=1), ncells, replace=FALSE)
	return(sce[,keep])
}
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

#' Run Gene Set Enrichment Analysis (GSEA) with MSigDb signatures
#'
#' This function performs GSEA on a named list of ranked genes, and restricts the analysis to specific MSigDb 
#' collections of gene signatures using the 'cats' and 'subcats' arguments.
#'
#' @param msigs A dataframe containing all MSigDb gene signatures.
#' @param ranked.genes A named list of ranked genes.
#' @param outdir The output directory for results.
#' @param outprefix The prefix for output files.
#' @param xlsx A list to store results that will be later written to an Excel file. Defaults to NULL.
#' @param cats A character vector specifying the main categories of gene sets to consider from MSigDb. Defaults to "H".
#' @param subcats A character vector specifying the subcategories of gene sets to consider from MSigDb. 
#'   Must match in length with 'cats'. Defaults to NULL.
#' @param ... Additional arguments to pass to the 'fgsea' function.
#'
#' @return A list of GSEA results if 'xlsx' is not NULL, otherwise, results are saved as files in the specified output directory.
#'
#' @examples
#' \dontrun{
#' runGSEA(msigs = msigdb, ranked.genes = my_genes, outdir = "./results", outprefix = "experiment1", 
#'         xlsx = list(), cats = c("H", "C3"), subcats = c("BP", "MIR"))
#' }
#' 
#' @note Ensure that 'cats' and 'subcats' vectors have equal lengths.
#' @note The function creates various output files including detailed GSEA results, 
#'   enrichment plots, and tables of top enriched pathways.
#' 
#' @author Jared Andrews
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

#' Run Gene Set Enrichment Analysis (GSEA) with custom signatures
#'
#' This function performs GSEA on a named list of ranked genes.
#'
#' @param sigs A named list of gene signatures.
#' @param ranked.genes A named list of ranked genes.
#' @param outdir The output directory for results.
#' @param outprefix The prefix for output files.
#' @param xlsx A list to store results that will be later written to an Excel file. Defaults to NULL.
#' @param ... Additional arguments to pass to the 'fgsea' function.
#'
#' @return A list of GSEA results if 'xlsx' is not NULL, otherwise, results are saved as files in the specified output directory.
#'
#' @examples
#' \dontrun{
#' runCustomGSEA(msigs = msigdb, ranked.genes = my_genes, outdir = "./results", outprefix = "experiment1", 
#'         xlsx = list())
#' }
#' 
#' @note The function creates various output files including detailed GSEA results, enrichment plots, and tables of top enriched pathways.
#' 
#' @author Jared Andrews
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
#' Summarize Gene Set Enrichment Analysis (GSEA) Results
#'
#' This function summarizes GSEA results by selecting the top significant gene sets. It creates an output directory
#' and saves the summarized results into this directory.
#'
#' @param gsea.list A list of GSEA results as returned by `runGSEA` or `runCustomGSEA`.
#' @param outdir The output directory for summarized results.
#' @param padj.th The significance threshold (adjusted p-value) for filtering gene sets. Defaults to 0.05.
#' @param top The number of top significant gene sets to consider. Defaults to 75.
#'
#' @return Invisible. The function creates an output directory and saves the summarized results there.
#'
#' @examples
#' \dontrun{
#' summarize_GSEA(gsea.list = my_gsea_results, outdir = "./summary", padj.th = 0.01, top = 50)
#' }
#'
#' @author Jared Andrews
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

#### GO Semantic Similarity Heatmaps

These cluster GO terms together by their similarity, allowing us to collapse closely related terms together. This is useful for summarizing the broad changes in each comparison.

```r
#' Simplify GSEA Results and Generate PDF
#'
#' This function simplifies GSEA (Gene Set Enrichment Analysis) results from
#' GO term genesets based on semantic similarity of the genesets.
#' It filters categorizes the results returned from fgsea based on the 
#' normalized enrichment score (NES) and adjust p-value.
#'
#' @param res Data frame containing the GSEA results from fgsea.
#' @param msig Data frame containing the gene set metadata.
#' @param outname Character string specifying the name of the output PDF file. Default is "simplified_GSEA.pdf".
#' @param pos.name Character string specifying the name for positively enriched sets. Default is "g1_enriched".
#' @param neg.name Character string specifying the name for negatively enriched sets. Default is "g2_enriched".
#' @param height Numeric specifying the height of the PDF. Default is 12.
#' @param width Numeric specifying the width of the PDF. Default is 12.
#' @param ... Additional parameters to pass to the `simplifyGOFromMultipleLists` function.
#'
#' @return NULL. The function generates a PDF file as a side effect.
#' @author Jared Andrews
#' 
#' @seealso \code{\link[simplifyEnrichment]{simplifyGOFromMultipleLists}}
simplify_GSEA <- function(res,
                          msig,
                          outname = "simplified_GSEA.pdf",
                          pos.name = "g1_enriched",
                          neg.name = "g2_enriched",
                          height = 12,
                          width = 12,
                          ...) {
  # Check for GO terms
  res$ID <- msig$gs_exact_source[match(res$path, msig$gs_name)]
  res$p.adjust <- res$padj
  
  if (nrow(res) < 5) {
	message("Not enough terms to cluster (<5), skipping.")
	return(NULL)
  }
  
  pos_res <- res[res$NES > 0, ]
  neg_res <- res[res$NES < 0, ]
  lt <- list()
  lt[[pos.name]] <- pos_res
  lt[[neg.name]] <- neg_res
  
  pdf(outname, height = height, width = width)
  simplifyGOFromMultipleLists(lt, ...)
  dev.off()
}

for (r in c("WT.midline.v.hemi.C5.GO.BP", "WT.midline.v.hemi.C5.GO.MF", "WT.midline.v.hemi.C5.GO.CC")) {
  fgsea.res <- xl.lists[[r]]
  simplify_GSEA(fgsea.res, msig, pos.name = "midline_enriched", neg.name = "hemispheric_enriched", 
			    outname = paste0("./GSEA/", r, ".simplifyHeatmap.pdf"), padj_cutoff = 0.05, 
			    min_term = 2, fontsize_range = c(7, 14), 
			    heatmap_param = list(col = c("blue", "white", "red"), breaks = c(1, 0.05, 0.0005)))
}
```
#### Plot Leading Edge Genes
GSEA returns the leading edge genes that are driving the score for a given signature. It can be useful to have a closer look at these genes in the form of boxplots and/or heatmaps.

```r
#' Plot Leading Edge Genes from GSEA
#'
#' This function plots the leading edge genes from a GSEA analysis, which are the core genes that contribute to 
#' the enrichment signal. These plots can help understand the gene expression patterns in the form of boxplots or heatmaps.
#'
#' @param dds A DESeqDataSet object.
#' @param gsea.lists A list of GSEA results as returned by `runGSEA` or `runCustomGSEA`.
#' @param annot.by A character string or vector for column name(s) in `colData(dds)` by which to annotate the samples
#' @param group.by A character string or vector for column name(s) in `colData(dds)` by which to group the samples
#' @param outdir The directory where the output plots should be saved.
#' @param use.assay A character string specifying the assay to use from the 'dds' object.
#' @param cells.use A character vector specifying the cells to include in the plot.
#' @param sig.thresh The significance threshold (adjusted p-value) for selecting gene sets. Defaults to 0.05.
#' @param group.by2 A secondary character string or vector for column name(s) in `colData(dds)` to further group the samples. Defaults to NULL.
#' @param split.by A character string or vector for column name(s) in `colData(dds)` to split the plot into multiple facets. Defaults to NULL.
#' @param swap.rownames A character string for `rowData` column to switch the rownames (e.g. "SYMBOL"). Defaults to NULL.
#' @param order.by A character string or vector for column name(s) in `colData(dds)` by which to order the samples Defaults to NULL.
#'
#' @return Invisible. The function saves the plots to the specified output directory.
#'
#' @examples
#' \dontrun{
#' plot_le(dds = my_dds, gsea.lists = my_gsea_results, annot.by = "group", group.by = "condition", 
#'         outdir = "./plots", use.assay = "counts", cells.use = c("cell1", "cell2"), 
#'         sig.thresh = 0.01, group.by2 = "timepoint", split.by = "treatment")
#' }
#'
#' @note The function may take a long time to execute if many gene sets are provided.
#'
#' @author Jared Andrews
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
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                          swap.rownames = swap.rownames)
          print(pl)
          pl <- dittoPlotVarsAcrossGroups(sce[,cells.use], le, group.by = group.by, 
                                          plots = c("vlnplot", "jitter", "boxplot"), 
                                          adjustment = "relative.to.max", assay = i, sub = i,
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                          swap.rownames = swap.rownames)
          print(pl)
          pl <- dittoPlotVarsAcrossGroups(sce[,cells.use], le, group.by = group.by, 
                                          plots = c("vlnplot", "jitter", "boxplot"), 
                                          adjustment = "none", assay = i, sub = i,
                                          vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                          swap.rownames = swap.rownames)
          print(pl)
          
          if (!is.null(group.by2) & !is.null(split.by)) {
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                            swap.rownames = swap.rownames)
            print(pl)
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), 
                                            adjustment = "relative.to.max", assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                            swap.rownames = swap.rownames)
            print(pl)
            pl <- dittoPlotVarsAcrossGroups(sce, le, group.by = group.by2, split.by = split.by,
                                            plots = c("vlnplot", "jitter", "boxplot"), 
                                            adjustment = "none", assay = i, sub = i,
                                            vlnplot.lineweight = 0.4, boxplot.lineweight = 0.5, 
                                            swap.rownames = swap.rownames)
            print(pl)
          }
        }
        
        dev.off()
        
        pdf(paste0(outdir, "/", ct, ".", path.name, ".heatmap.pdf"), width = 5, height = 7)
        for (i in use.assay) {
          pl <- dittoHeatmap(sce, le, annot.by = annot.by, cells.use = cells.use, show_colnames = FALSE,
                             breaks = seq(-3, 3, length.out = 51), cluster_rows = FALSE,
                             fontsize_row = 6, cluster_cols = FALSE, assay = i, sub = i, 
                             swap.rownames = swap.rownames)
          grid.draw(pl)
          
          pl <- dittoHeatmap(sce, le, annot.by = annot.by, show_colnames = FALSE,
                             breaks = seq(-3, 3, length.out = 51), cluster_rows = FALSE,
                             fontsize_row = 6, cluster_cols = FALSE, assay = i, sub = i, 
                             swap.rownames = swap.rownames)
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
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("stringi")
library("pathview")
library("ReactomePA")

#' @param res.list Named list of DESeq2 results data.frames.
#' @param padj.th Numeric scalar to use as significance threshold for DE genes.
#' @param lfc.th Numeric scalar to use as log fold change threshold for DE genes.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.col Character scalar indicating name of gene ID column for each data.frame in \code{res.list}
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param organism Character scalar indicating species in KEGG format ("hsa", "mmu", etc).
#' @param ... Passed to \code{compareCluster}.
#' @author Jared Andrews
run_enrichKEGG <- function(res.list, padj.th = 0.05, lfc.th = 0, outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.col = "ENSEMBL", id.type = "ENSEMBL", 
                         organism = "hsa", ...) {
  # Do GO enrichment on up/downregulated genes.
  for (r in names(res.list)) {
    df <- res.list[[r]]

	if (lfc.th != 0) {
      r <- paste0(r,"-LFC", round(lfc.th, digits = 3), "filt")
    }

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
                  downregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange < -lfc.th],
                  all_de = df[[id.col]][df$padj < padj.th])
    
    skip <- FALSE
									
    tryCatch({
      genes$upregulated <- bitr(genes$upregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$downregulated <- bitr(genes$downregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$all_de <- bitr(genes$all_de, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
    }, 
    error = function(e) {
      message("There was an error: ", e)
      message("Most likely, no gene identifiers for hits could be mapped to entrez IDs.")
      message("Proceeding to next comparison in results list.")
      skip <<- TRUE
    })

	if (skip) {
      next
    }
    
    # Remove lowly expressed genes.
    bg <- df[[id.col]][!is.na(df$padj)]
    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    
    ck <- compareCluster(geneCluster = genes, fun = enrichKEGG, keyType = "kegg", universe = bg, organism = organism, ...)
    if (!is.null(ck)) {
      ck <- setReadable(ck, OrgDb = OrgDb, keyType="ENTREZID")
      
      # Term similarities via Jaccard Similarity index.
      ego <- pairwise_termsim(ck)
      
      height = 4 + (0.015 * length(ego@compareClusterResult$Cluster))
      
      pdf(paste0(out, "/KEGG_Enrichments.Top20_perGroup.pdf"), width = 6, height = height)
      p <- dotplot(ego, showCategory = 20, font.size = 7)
      print(p)
      p <- dotplot(ego, size = "count", showCategory = 20, font.size = 7)
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
                   species    = organism,
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

run_enrichKEGG(res, OrgDb = "org.Mm.eg.db", organism = "mmu")
```

#### Reactome Enrichment

```r
#' @param res.list Named list of DESeq2 results data.frames.
#' @param padj.th Numeric scalar to use as significance threshold for DE genes.
#' @param lfc.th Numeric scalar to use as log fold change threshold for DE genes.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.col Character scalar indicating name of gene ID column for each data.frame in \code{res.list}
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param organism Character scalar indicating ReactomePA-supported species ("human", "mouse").
#' @param ... Passed to \code{compareCluster}.
#' @author Jared Andrews
run_enrichPathway <- function(res.list, padj.th = 0.05, lfc.th = 0, outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.col = "ENSEMBL", id.type = "ENSEMBL", organism = "human", ...) {
  # Do GO enrichment on up/downregulated genes.
  for (r in names(res.list)) {
    df <- res.list[[r]]

	if (lfc.th != 0) {
      r <- paste0(r,"-LFC", round(lfc.th, digits = 3), "filt")
    }

    out <- file.path(outdir, r)
    dir.create(paste0(out, "/Reactome_pathways"), showWarnings = FALSE, recursive = TRUE)
    
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
    
    genes <- list(upregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange > lfc.th],
                  downregulated = df[[id.col]][df$padj < padj.th  & df$log2FoldChange < -lfc.th],
                  all_de = df[[id.col]][df$padj < padj.th])
    
    skip <- FALSE
									
    tryCatch({
      genes$upregulated <- bitr(genes$upregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$downregulated <- bitr(genes$downregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$all_de <- bitr(genes$all_de, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
    }, 
    error = function(e) {
      message("There was an error: ", e)
      message("Most likely, no gene identifiers for hits could be mapped to entrez IDs.")
      message("Proceeding to next comparison in results list.")
      skip <<- TRUE
    })

	if (skip) {
      next
    }
    
    # Remove duplicate IDs.
    gl <- gl[unique(names(gl))]
    
    # Remove lowly expressed genes.
    bg <- df[[id.col]][!is.na(df$padj)]
    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    
    ck <- compareCluster(geneCluster = genes, fun = enrichPathway, universe = bg, 
                         readable = TRUE, organism = organism, ...)
    if (!is.null(ck)) {
      # Term similarities via Jaccard Similarity index.
      ego <- pairwise_termsim(ck)
      
      # Adjust plot height based on number of terms.
      height = 4 + (0.015 * length(ego@compareClusterResult$Cluster))
      
      pdf(paste0(out, "/Reactome_Enrichments.Top20_perGroup.pdf"), width = 6, height = height)
      p <- dotplot(ego, showCategory = 20, font.size = 7)
      print(p)
      p <- dotplot(ego, size = "count", showCategory = 20, font.size = 7)
      print(p)
      dev.off()
      
      pdf(paste0(out, "/Reactome_Enrichments.termsim.Top20_perGroup.pdf"), width = 9, height = 9)
      p <- emapplot(ego, pie="count", cex_category=0.9, cex_label_category = 0.9, 
                    layout="kk", repel = TRUE, showCategory = 20)
      print(p)
      dev.off()
      
      if (nrow(ego) > 2) {
        pdf(paste0(out, "/Reactome_Enrichments.termsim.Top30_Tree.pdf"), width = 17, height = 14)
        p <- treeplot(ego, showCategory = 30, fontsize = 4, offset.params = list(bar_tree = rel(2.5), tiplab = rel(2.5), extend = 0.3, hexpand = 0.1), cluster.params = list(method = "ward.D", n = min(c(6, ceiling(sqrt(nrow(ego))))), color = NULL, label_words_n = 5, label_format = 30))
        print(p)
        dev.off()
        
        pdf(paste0(out, "/Reactome_Enrichments.termsim.Top10_FullNet.pdf"), width = 15, height = 15)
        p <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
        print(p)
        dev.off()
        
        pdf(paste0(out, "/Reactome_Enrichments.termsim.Top5_FullNet.pdf"), width = 12, height = 12)
        p <- cnetplot(ego, showCategory = 5, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
        print(p)
        dev.off()
      }
      
      # Network plot for each pathway with FC values.
      # for (x in ego@compareClusterResult$Description) {
      #   # Some pathways have backslashes, which will break file creation.
      #   x_out <- str_replace_all(x, "/", "_")
      #   
      #   # For when it inevitably wants to crash due to not finding a pathway name or such.
      #   tryCatch(
      #     {
      #       pdf(paste0(out, "/Reactome_pathways/", x_out, ".pdf"), width = 11, height = 11)
      #       p <- viewPathway(x, readable = TRUE, foldChange = gl, organism = organism)
      #       vals <- p$data$color[!is.na(p$data$color)]
      #       l <- max(abs(as.numeric(vals)))
      #       p <- p + scale_color_gradient2(limits = c(-l,l), mid = "grey90", 
      #                                      high = "red", low = "navyblue")
      #       print(p)
      #       dev.off()
      #       dev.off()
      #       dev.off()
      #       dev.off()
      #       dev.off()
      #       dev.off()
      #     },
      #     error = function(e) {"bah"}
      #   )
      # }
      
      saveRDS(ego, file = paste0(out, "/enrichPathway.reactome.RDS"))
      ego <- as.data.frame(ego)
      write.table(ego, file = paste0(out, "/enrichPathway.reactome.results.txt"), 
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

run_enrichPathway(res, OrgDb = "org.Mm.eg.db", organism = "mouse")
```

#### GO Enrichment

```r
#' @param res.list Named list of DESeq2 results data.frames.
#' @param padj.th Numeric scalar to use as significance threshold for DE genes.
#' @param lfc.th Numeric scalar to use as log fold change threshold for DE genes.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.col Character scalar indicating name of gene ID column for each data.frame in \code{res.list}
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param onts Character vector indicating ontologies to test individually. 
#'   Options must be one or more of "ALL", "BP", "CC", or "MF". 
#'   Default uses all of those options. 
#' @param ... Passed to \code{compareCluster}
#' @author Jared Andrews
run_enrichGO <- function(res.list, padj.th = 0.05, lfc.th = 0, outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.col = "ENSEMBL", id.type = "ENSEMBL", 
                         onts = c("BP", "MF", "CC", "ALL"), ...) {
  # Do GO enrichment on up/downregulated genes.
  for (r in names(res.list)) {
    df <- res.list[[r]]

	if (lfc.th != 0) {
      r <- paste0(r,"-LFC", round(lfc.th, digits = 3), "filt")
    }

    out <- file.path(outdir, r)
    dir.create(paste0(out, "/GO_enrichments"), showWarnings = FALSE, recursive = TRUE)
    
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
    gl <- sort(gl, decreasing = TRUE)
    # Remove duplicate IDs.
    gl <- gl[unique(names(gl))]
    
    genes <- list(upregulated = df[[id.col]][df$padj < padj.th & df$log2FoldChange > lfc.th],
                  downregulated = df[[id.col]][df$padj < padj.th  & df$log2FoldChange < -lfc.th],
                  all_de = df[[id.col]][df$padj < padj.th])

	skip <- FALSE
									
    tryCatch({
      genes$upregulated <- bitr(genes$upregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$downregulated <- bitr(genes$downregulated, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
  
      genes$all_de <- bitr(genes$all_de, fromType = id.type, toType = "ENTREZID", 
                            OrgDb = OrgDb)$ENTREZID
    }, 
    error = function(e) {
      message("There was an error: ", e)
      message("Most likely, no gene identifiers for hits could be mapped to entrez IDs.")
      message("Proceeding to next comparison in results list.")
      skip <<- TRUE
    })

	if (skip) {
      next
    }
    
    # Remove lowly expressed genes.
    bg <- df[[id.col]][!is.na(df$padj)]
    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    
    for (ont in onts) {
    
      ck <- compareCluster(geneCluster = genes, fun = enrichGO, universe = bg, 
                           readable = TRUE, ont = ont, OrgDb = OrgDb, ...)
      
      if (!is.null(ck)) {
        # Term similarities via Jaccard Similarity index.
        ego <- pairwise_termsim(ck)
        
        # Adjust plot height based on number of terms.
        height = 3.75 + (0.025 * length(ego@compareClusterResult$Cluster))
        
        pdf(paste0(out, "/GO_Enrichments.Top20_perGroup.", ont, ".pdf"), width = 6, height = height)
        p <- dotplot(ego, showCategory = 20, font.size = 7)
        print(p)
        p <- dotplot(ego, size = "count", showCategory = 20, font.size = 7)
        print(p)
        dev.off()
        
        pdf(paste0(out, "/GO_Enrichments.termsim.Top20_perGroup.", ont, ".pdf"), width = 9, height = 9)
        p <- emapplot(ego, pie="count", cex_category=0.9, cex_label_category = 0.9,
                      layout="kk", repel = TRUE, showCategory = 20)
        print(p)
        dev.off()
        
        
        if (nrow(ego) > 2) {
          pdf(paste0(out, "/GO_Enrichments.termsim.Top30_Tree.", ont, ".pdf"), width = 17, height = 14)
          p <- treeplot(ego, showCategory = 30, fontsize = 4, offset.params = list(bar_tree = rel(2.5), tiplab = rel(2.5), extend = 0.3, hexpand = 0.1), cluster.params = list(method = "ward.D", n = min(c(6, ceiling(sqrt(nrow(ego))))), color = NULL, label_words_n = 5, label_format = 30))
          print(p)
          dev.off()
          
          pdf(paste0(out, "/GO_Enrichments.termsim.Top10_FullNet.", ont, ".pdf"), width = 15, height = 15)
          p <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
          print(p)
          dev.off()
          
          pdf(paste0(out, "/GO_Enrichments.termsim.Top5_FullNet.", ont, ".pdf"), width = 13, height = 13)
          p <- cnetplot(ego, showCategory = 5, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
          print(p)
          dev.off()
        }
        
        saveRDS(ego, file = paste0(out, "/enrichGO.", ont, ".RDS"))
        ego <- as.data.frame(ego)
        write.table(ego, file = paste0(out, "/enrichGO.results.", ont, ".txt"), 
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  }
}

run_enrichGO(res, OrgDb = "org.Mm.eg.db")
```

#### KEGG/Reactome Enrichment (simple)

This version just uses lists of genes that can be defined however which way.

```r
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("stringi")
library("pathview")
library("ReactomePA")

#' @param genes Character vector of gene IDs to test for enrichment.
#' @param bg Character vector of gene IDs to be used as background.
#' @param name Character scalar for output name.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param kegg.organism Character scalar indicating species in KEGG format ("hsa", "mmu", etc).
#' @param reactome.organism Character scalar indicating species in Reactome format ("human", "mouse", etc).
#' @param fun Character scalar indicating which "enrich" function to run ("enrichKEGG", "enrichPathway").
#' @param ... Passed to specified enrichments function.
#' @author Jared Andrews
run_enrich_simple <- function(genes, bg, name = "sample", outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.type = "ENSEMBL", 
                         kegg.organism = "hsa", reactome.organism = "human", fun = "enrichKEGG", ...) {

  out <- file.path(outdir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  
  # Strip gene version info if ensembl.
  if (id.type == "ENSEMBL") {
    xx <- strsplit(genes, "\\.")
    genes <- unlist(lapply(xx, FUN = function(x) x[1]))
    
    xx <- strsplit(bg, "\\.")
    bg <- unlist(lapply(xx, FUN = function(x) x[1]))
  }

  skip <- FALSE
						
  tryCatch({
    genes <- bitr(genes, fromType = id.type, toType = "ENTREZID", 
                          OrgDb = OrgDb)$ENTREZID
  }, 
  error = function(e) {
    message("There was an error: ", e)
    message("Most likely, no gene identifiers for hits could be mapped to entrez IDs.")
    skip <<- TRUE
  })
  
  if (skip) {
	next
  }
  
  bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
  
  if (fun == "enrichKEGG") {
    ego <- enrichKEGG(gene = genes, universe = bg, organism = kegg.organism, ...)
  } else if (fun == "enrichPathway") {
    ego <- enrichPathway(gene = genes, universe = bg, organism = reactome.organism, ...)
  } 

  if (nrow(as.data.frame(ego)) > 0) {
    # Term similarities via Jaccard Similarity index.
    ego <- pairwise_termsim(ego)
    ego <- setReadable(ego, OrgDb = OrgDb, keyType="ENTREZID")
    pdf(paste0(out, "/", fun, ".Top30.pdf"), width = 6, height = 8)
    p <- dotplot(ego, showCategory = 30, font.size = 7)
    print(p)
    p <- barplot(ego, showCategory = 30, font.size = 7)
    print(p)
    dev.off()
    
    if (nrow(as.data.frame(ego)) > 2) {
      pdf(paste0(out, "/", fun, ".termsim.Top30_Tree.pdf"), width = 17, height = 14)
      p <- treeplot(ego, showCategory = 30, fontsize = 4, offset.params = list(bar_tree = rel(2.5), tiplab = rel(2.5), extend = 0.3, hexpand = 0.1), cluster.params = list(method = "ward.D", n = min(c(6, ceiling(sqrt(nrow(ego))))), color = NULL, label_words_n = 5, label_format = 30))
      print(p)
      dev.off()
      
      pdf(paste0(out, "/", fun, ".termsim.Top10_FullNet.pdf"), width = 15, height = 15)
      p <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
      print(p)
      dev.off()
      
      pdf(paste0(out, "/", fun, ".termsim.Top5_FullNet.pdf"), width = 12, height = 12)
      p <- cnetplot(ego, showCategory = 5, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
      print(p)
      dev.off()
    }
    
    saveRDS(ego, file = paste0(out, "/", fun, ".results.RDS"))
    ego <- as.data.frame(ego)
    write.table(ego, file = paste0(out, "/", fun, ".results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
```

#### GO Enrichment (simple)

This version just uses lists of genes that can be defined however which way.

```r
#' @param genes Character vector of gene IDs to test for enrichment.
#' @param bg Character vector of gene IDs to be used as background.
#' @param name Character scalar for output name.
#' @param outdir Character scalar for output directory.
#' @param OrgDb Character scalar for annotation database to use.
#' @param id.type Character scalar indicating type of gene ID used. See \code{keytypes(org.Hs.eg.db)} for all options.
#' @param onts Character vector indicating ontologies to test individually. 
#'   Options must be one or more of "ALL", "BP", "CC", or "MF". 
#'   Default uses all of those options. 
#' @param ... Passed to specified enrichments function.
#' @author Jared Andrews
run_enrichGO_simple <- function(genes, bg, name = "sample", outdir = "./enrichments",
                         OrgDb = "org.Hs.eg.db", id.type = "ENSEMBL", onts = c("BP", "MF", "CC", "ALL"), ...) {

  out <- file.path(outdir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  
  # Strip gene version info if ensembl.
  if (id.type == "ENSEMBL") {
    xx <- strsplit(genes, "\\.")
    genes <- unlist(lapply(xx, FUN = function(x) x[1]))
    
    xx <- strsplit(bg, "\\.")
    bg <- unlist(lapply(xx, FUN = function(x) x[1]))
  }

  skip <- FALSE
						
  tryCatch({
    genes <- bitr(genes, fromType = id.type, toType = "ENTREZID", 
                          OrgDb = OrgDb)$ENTREZID

  }, 
  error = function(e) {
    message("There was an error: ", e)
    message("Most likely, no gene identifiers for hits could be mapped to entrez IDs.")
    skip <<- TRUE
  })

  if (skip) {
	next
  }
  
  bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID

  for (ont in onts) {
    ego <- enrichGO(genes, OrgDb = OrgDb, universe = bg, ont = ont, readable = TRUE, ...)

    if (nrow(as.data.frame(ego)) > 0) {
      # Term similarities via Jaccard Similarity index.
      ego <- pairwise_termsim(ego)
      pdf(paste0(out, "/enrichGO.", ont, ".Top30.pdf"), width = 6, height = 8)
      p <- dotplot(ego, showCategory = 30, font.size = 7)
      print(p)
      p <- barplot(ego, showCategory = 30, font.size = 7)
      print(p)
      dev.off()
      
      if (nrow(as.data.frame(ego)) > 2) {
        pdf(paste0(out, "/enrichGO.", ont, ".termsim.Top30_Tree.pdf"), width = 17, height = 14)
        p <- treeplot(ego, showCategory = 30, fontsize = 4, offset.params = list(bar_tree = rel(2.5), tiplab = rel(2.5), extend = 0.3, hexpand = 0.1), cluster.params = list(method = "ward.D", n = min(c(6, ceiling(sqrt(nrow(ego))))), color = NULL, label_words_n = 5, label_format = 30))
        print(p)
        dev.off()
        
        pdf(paste0(out, "/enrichGO.", ont, ".termsim.Top10_FullNet.pdf"), width = 15, height = 15)
        p <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
        print(p)
        dev.off()
        
        pdf(paste0(out, "/enrichGO.", ont, ".termsim.Top5_FullNet.pdf"), width = 12, height = 12)
        p <- cnetplot(ego, showCategory = 5, cex.params = list(category_label = 1.3, gene_label = 0.9, category_node = 1, gene_node = 1), layout = "kk")
        print(p)
        dev.off()
      }
      
      saveRDS(ego, file = paste0(out, "/enrichGO.", ont, ".results.RDS"))
      ego <- as.data.frame(ego)
      write.table(ego, file = paste0(out, "/enrichGO.", ont, ".results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

rezzies <- c("K_0.v.W_0-b_cul", "K_1.v.W_1-b_cul", "K_2.v.W_2-b_cul", "K_5.v.W_5-b_cul", "K_7.v.W_7-b_cul")
bgs <- list(bg0 = res$`K_0.v.W_0-b_cul`$SYMBOL[!is.na(res$`K_0.v.W_0-b_cul`$padj)],
  bg1 = res$`K_1.v.W_1-b_cul`$SYMBOL[!is.na(res$`K_1.v.W_1-b_cul`$padj)],
  bg2 = res$`K_2.v.W_2-b_cul`$SYMBOL[!is.na(res$`K_2.v.W_2-b_cul`$padj)],
  bg5 = res$`K_5.v.W_5-b_cul`$SYMBOL[!is.na(res$`K_5.v.W_5-b_cul`$padj)],
  bg7 = res$`K_7.v.W_7-b_cul`$SYMBOL[!is.na(res$`K_7.v.W_7-b_cul`$padj)])

for (i in seq_along(df_lists)) {
  n <- names(df_lists)[i]
  bg <- ifelse(grepl("0d", n), "bg0", 
               ifelse(grepl("1d", n), "bg1",
               ifelse(grepl("2d", n), "bg2",
               ifelse(grepl("5d", n), "bg5",
               ifelse(grepl("7d", n), "bg7")))))
  g <- df_lists[[i]]
  
  run_enrich_simple(g, bg = bgs[[bg]], OrgDb = "org.Mm.eg.db", kegg.organism = "mmu", reactome.organism = "mouse", id.type = "SYMBOL", outdir = "./eed_comparison/enrichments", name = n, fun = "enrichKEGG")
  
  run_enrich_simple(g, bg = bgs[[bg]], OrgDb = "org.Mm.eg.db", kegg.organism = "mmu", reactome.organism = "mouse", id.type = "SYMBOL", outdir = "./eed_comparison/enrichments", name = n, fun = "enrichPathway")
  
  run_enrichGO_simple(g, bg = bgs[[bg]], OrgDb = "org.Mm.eg.db", id.type = "SYMBOL", outdir = "./eed_comparison/enrichments", name = n)
}
```

### CNV Calling from Methylation Array
This spits out typical genome-wide CNV plots, segmentation files, bins, and IGV tracks from Illumina methylation arrays. Users can add details regions for labels if they'd like. When mixing both 450k and EPIC arrays, set `array_type = "overlap"`.
```r
library("minfi")
library("conumee")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

#' @param meta Character scalar for the path to samplesheet with sample metadata, each row containing a sample.
#' @param basedir Character scalar for the base directory containing IDATs.
#' @param controls Character vector for control sample names found in \code{name_col}.
#' @param outdir Character scalar for output directory.
#' @param chr Character vector indicating chromosomes to plot in genome plots. If provided, an additional
#'   set of genome plots will be created with only these chromosomes.
#' @param exclude_regions GRanges object containing regions to exclude from the CN plots.
#' @param detail_regions GRanges object containing regions to label.
#' @param array_type Character scalar indicating more array type. Options are "450k", "EPIC", or "overlap" for
#'   datasets with mixed arrays.
#' @param idat_cols Character scalar or vector containing column names for sample identifier columns to paste together.
#'   This should be the IDAT ID.
#' @param name_col Character scalar for column containing sample names.
run_conumee_CNV <- function(meta, basedir, controls, outdir, chr = "all", exclude_regions = NULL, detail_regions = NULL, 
                            array_type = "450k", idat_cols = c("Sentrix_ID", "Sentrix_Position"),
                            name_col = "Sample") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ## Load data.
  meta <- read.csv(meta)
  meta$Basename <- file.path(basedir, apply(meta[,idat_cols, drop = FALSE], MARGIN = 1, FUN = paste0, collapse = "_"))

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
    x <- CNV.segment(x)
    
    CNV.genomeplot(x, main = s)
    
    x <- CNV.detail(x)
    CNV.genomeplot(x, main = paste0(s, " - detailed"))
    
    if (length(chr) > 1) {
      CNV.genomeplot(x, main = paste0(s, " - ", paste(chr, collapse = ", ")), chr = chr)
    }
    
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

### Find Common Elements from Multiple Vectors
```r
Reduce(intersect, list(a,b,c))
```



### Differential Gene Expression via `DESeq2`
This is a super lazy function to run through a list of contrasts and create differential gene expression results for each.

```r
#' Get DESeq2 Results
#'
#' This function obtains a set of comparisons from a DESeq2 analysis, given a named list of contrasts. It allows additional model 
#' parameters to be specified and a design matrix to be manually adjusted. 
#'
#' @param dds An object of class DESeqDataSet.
#' @param res.list A list of DESeq2 result tables. Allows the fuction to be run multiple times if needed and append to the same list.
#' @param contrasts A named list of contrasts.
#' @param user.mat A logical indicating whether a user-specified model matrix is provided. Defaults to FALSE.
#' @param block A vector of additional terms to be considered in the model, beyond the main effect. Defaults to NULL.
#' @param design The design formula or matrix. If a matrix is provided, ensure 'user.mat' is set to TRUE. Defaults to NULL.
#' @param alpha The significance level for hypothesis testing. Defaults to 0.05.
#' @param lfc.th A numeric vector of log2 fold-change thresholds. Defaults to c(log2(1.5), log2(2)).
#' @param shrink.method The method used for shrinkage estimation. Defaults to "apeglm".
#' @param outdir The directory where the output should be saved. Defaults to "./de".
#' @param BPPARAM The BiocParallelParam object specifying the parallel back-end to be used. Defaults to NULL.
#' 
#' @return A list of DESeq2 result tables for the specified contrasts, saved to the specified output directory.
#' 
#' @examples
#' \dontrun{
#' get_DESEQ2_res(dds, res.list, contrasts, user.mat = TRUE, block = c("term1", "term2"), 
#'                design = my_design, alpha = 0.01, lfc.th = c(log2(2), log2(3)), 
#'                shrink.method = "normal", outdir = "./my_results", BPPARAM = MulticoreParam(2))
#' }
get_DESEQ2_res <- function(dds, res.list, contrasts, user.mat = FALSE, block = NULL, design = NULL, alpha = 0.05, 
						   lfc.th = c(log2(1.5), log2(2)), shrink.method = "apeglm", outdir = "./de", BPPARAM = NULL) {
  
  dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)
  
  for (i in seq_along(contrasts)) {
    rname <- names(contrasts)[i]
    
    # If user-supplied matrix, contrast must be in list format.
    if (user.mat) {
      con <- contrasts[[i]]
      message("Setting shrink.method to 'ashr' to work with list contrasts due to user-specified model matrix.")
      shrink.method <- "ashr"
    } else {
      con <- contrasts[[i]]
      coef <- paste(con[1], con[2], "vs", con[3], sep = "_")
      dds[[con[1]]] <- relevel(dds[[con[1]]], ref = con[3])
    }
    
    if (!is.null(design)) {
      desgn <- design
    } else if (!is.null(block)) {
      desgn <- as.formula(paste0("~",paste0(c(block, con[1]), collapse = "+")))
    } else {
      desgn <- as.formula(paste0("~", con[1]))
    }
    
    message(paste0("Design for ", paste(con[1], con[2], "vs", con[3], sep = "_"),
                   " is ", paste0(as.character(desgn))))
    
    dds <- DESeqDataSet(dds, design = desgn)
    dds <- DESeq(dds, BPPARAM = BPPARAM)

    res1 <- results(dds, contrast = con, alpha = alpha)
    res1$ENSEMBL <- rownames(res1)
    res1$SYMBOL <- rowData(dds)$SYMBOL
    
    if (!is.null(shrink.method)) {
      out.name <- paste0(rname, "-shLFC")
      
      # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
      if (shrink.method == "ashr") {
        coef <- NULL
      }
      shrink <- lfcShrink(dds, res = res1, coef = coef, type = shrink.method)
      shrink$ENSEMBL <- rownames(shrink)
      shrink$SYMBOL <- rowData(dds)$SYMBOL
      rownames(shrink) <- shrink$SYMBOL
      shrink <- as.data.frame(shrink)
      res.list[[out.name]] <- shrink
      write.table(shrink, file = paste0(outdir, "/", rname, ".shrinkFC.padj.", alpha, ".txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    
    rownames(res1) <- res1$SYMBOL
    res1 <- as.data.frame(res1)
    res.list[[rname]] <- res1
    
    write.table(res1, file = paste0(outdir, "/", rname, ".padj.", alpha, ".txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    
    for (l in lfc.th) {
      
      res <- results(dds, contrast = con, alpha = alpha, lfcThreshold = l)
      res$ENSEMBL <- rownames(res)
      res$SYMBOL <- rowData(dds)$SYMBOL
      
      if (!is.null(shrink.method)) {
        # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
        if (shrink.method == "ashr") {
          coef <- NULL
        }
        out.name <- paste0(rname, "-shLFC", l)
        shrink <- lfcShrink(dds, res = res, coef = coef, type = shrink.method)
        shrink$ENSEMBL <- rownames(shrink)
        shrink$SYMBOL <- rowData(dds)$SYMBOL
        rownames(shrink) <- shrink$SYMBOL
        shrink <- as.data.frame(shrink)
        res.list[[out.name]] <- shrink
        write.table(shrink, file = paste0(outdir, "/", rname, ".shrinkLFC_thresh.", l,".padj.", alpha, ".txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      }
      
      rownames(res) <- res$SYMBOL
      out.name <- paste0(rname, "-LFC", l)
      res <- as.data.frame(res)
      res.list[[out.name]] <- res
      
      write.table(res, file = paste0(outdir, "/", rname, ".LFC_thresh.", l,".padj.", alpha, ".txt"), 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
  
  return(res.list)
}

res <- list()

contrasts = list("shPDGFRA_2.v.shScr" = c("shRNA", "shPDGFRA_2", "shScr"),
                 "shPDGFRA_3.v.shScr" = c("shRNA", "shPDGFRA_3", "shScr"),
                 "shZFP36L1_1.v.shScr" = c("shRNA", "shZFP36L1_1", "shScr"),
                 "shZFP36L1_2.v.shScr" = c("shRNA", "shZFP36L1_2", "shScr"),
                 "shPTPRZ1_3.v.shScr" = c("shRNA", "shPTPRZ1_3", "shScr"))
                 
res <- get_DESEQ2_res(dds, res.list = res, contrasts = contrasts)
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

#### Concatenate FASTQs for Given Sample

For merging reads across multiple lanes, etc.

```bash
function concat_fastq {

    sample=$1

    if [[ -z "$sample" ]]; then
        echo "A sample name must be provided as an argument."
        return 1
    fi


    # Check if files for R2 exist
    if ls "${sample}"_L00*_R2_001.fastq.gz 1> /dev/null 2>&1; then
        # This is a paired-end sample
        echo "Processing paired-end sample: $sample"

        # Concatenate R1
        cat "${sample}"_L00*_R1_001.fastq.gz > "${sample}_merged_R1_001.fastq.gz"
        if [[ $? -ne 0 ]]; then
            echo "An error occurred while concatenating R1 files."
            return 1
        fi

        # Concatenate R2
        cat "${sample}"_L00*_R2_001.fastq.gz > "${sample}_merged_R2_001.fastq.gz"

        if [[ $? -ne 0 ]]; then
            echo "An error occurred while concatenating R2 files."
            return 1
        fi
    else

        # This is a single-end sample
        echo "Processing single-end sample: $sample"

        # Concatenate R1 only
        cat "${sample}"_L00*_R1_001.fastq.gz > "${sample}_merged_R1_001.fastq.gz"
        if [[ $? -ne 0 ]]; then
            echo "An error occurred while concatenating R1 files."
            return 1
        fi
    fi

    echo "Concatenation complete for sample: $sample"
}
```

To use from a directory of FASTQs:
```bash
for sample_name in $(ls *_L00*_R1_001.fastq.gz | rev | cut -d'_' -f4- | rev | sort | uniq); do
    concat_fastq "$sample_name"
done
```

#### Recursively Get File Paths with Suffix

Useful for grabbing FASTQs from seq runs and all.

```bash
# Function to recursively get all files with a given suffix from a directory, 
# sort them, and output their full paths to a text file. 
# Usage: get_files_with_suffix_sorted <directory> <suffix> <output_file> 
function get_files_with_suffix { 
	local directory=$1 
	local suffix=$2 
	local output_file=$3 
	
	# Ensure the directory exists 
	if [[ ! -d "$directory" ]]; then 
		echo "The specified directory does not exist: $directory" 
		return 1 
	fi 
	
	# Use find to get files with the given suffix, sort them, and output to the specified file 
	find "$directory" -type f -name "*.$suffix" | sort > "$output_file" 
} 

# Example usage: # get_files_with_suffix /path/to/directory txt output_sorted.txt
```

### `tar` and `gzip` Directory

```bash
tar czf name_of_archive_file.tar.gz name_of_directory_to_tar
```

### Remove All Characters Up to Delimiter
Includes first instance of delimiter. ':' is the delimiter in this case.

```bash
sed 's/^[^:]*://' file
```

Useful for file renaming as well, e.g.:
```bash
for f in *.idat.gz; do samp=$(echo "$f" | sed 's/^[^_]*_//'); mv "$f" "$samp"; done
```

### Kill All LSF Jobs with String in Name

For when you screw up an array or whatnot.

```bash
 bjobs -w | grep 'bbsplit50' | awk '{print $1}' | xargs bkill
```

## File Manipulation Tasks

> If you don't love data cleaning, you don't love bioinformatics.
> -some guy

### Randomly Downsample BAM to Set Number of Reads/Read Pairs

For instance, to 5 million here. It will never keep a read but not its mate. From [biostars](https://www.biostars.org/p/9485840/#9485906).

```bash
reads=5000000
bam=your.bam
fraction=$(samtools idxstats $bam | cut -f3 | awk -v ct=$reads 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view -b -s ${fraction} foo.bam > sampled.bam
```

### BAM to FASTQ
Gimme them reads back. Note that these may not be identical to the original FASTQ files if unaligned reads weren't retained in the BAM file. Note this is for an IBM LSF cluster, just yoink the `bamtofastq` command if you're running locally.

```bash
#!/bin/bash
module load biobambam

for f in *.bam; do

    base="${f%%.bam}"

    bsub -P xeno -J bambam -B -N -n 4 -R "rusage[mem=4GB]" -M 4GB -q priority "bamtofastq F=$base.1.fastq.gz F2=$base.2.fastq.gz \
	    gz=1 filename=$f;"

done
```

### Extract Random Reads from FASTQ

For checking or whatnot. Change `samplerate` to extract a given percentage of reads, set to 0.1% here. `#` is used by bbmap to fill in 1 & 2, so this is for a paired end sample and will return two output files. Can also use `in1`, `in2`, `out1`, and `out2` if you'd prefer.

```bash
module load bbmap
reformat.sh in=thefastq_R#_001.fastq.gz out=thefastq_R#_001.subsampled.fastq.gz samplerate=0.001

# Or if a specific number of reads (pairs) are desired, e.g. 61,111,111 for 5x coverage with paired 150 bp read WGS.
reformat.sh in=thefastq_R#_001.fastq.gz out=thefastq_R#_001.subsampled.fastq.gz samplereadstarget=61111111
# 10x coverage with 150 bp reads.
reformat.sh in=thefastq_R#_001.fastq.gz out=thefastq_R#_001.subsampled.fastq.gz samplereadstarget=122222222

# 5x coverage with paired 100 bp read WGS.
reformat.sh in=thefastq_R#_001.fastq.gz out=thefastq_R#_001.subsampled.fastq.gz samplereadstarget=91666666
# 10x coverage with 100 bp reads.
reformat.sh in=thefastq_R#_001.fastq.gz out=thefastq_R#_001.subsampled.fastq.gz samplereadstarget=183333333
```

### Remove N Characters from End of Field in CSV

For removing PAM sequence, etc.

```bash
awk -F, -v OFS=',' '{ $2=substr($2, 1, length($2)-3) } 1' input.csv > output.csv
```