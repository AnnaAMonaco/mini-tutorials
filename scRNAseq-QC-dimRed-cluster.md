# WORK IN PROGRESS
Here you will learn the basics of analysing scRNA-seq datasets in R. I will cover the following steps:
- Loading mapped data from salmon alevin and from cellranger
- Quality control and filtering
- Integrating multiple datasets
- Normalisation
- Dimensionality reduction and clustering
- Marker genes and cell type annotation

These are the libraries you will need to have installed and loaded:
```{r}
# single-cell related
library(fishpond)
library(tximport)
library(tximeta)
library(SingleCellExperiment)
library(Seurat)
library(scran)
library(scater)
library(scDblFinder)
library(BiocSingular)
# data wrangling helpers
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
# plotting helpers
library(cowplot)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(viridis)
```

## Load data
If you quantified with `alevin-fry` you can load the dataset with `loadFry()`. The resulting object is a `SingleCellExperiment` object, so it will also need to be transformed into a `SeuratObject` for downstream processing.
```{r}
fryDir <- read.delim("samples.txt", header=FALSE)
fryDir <- as.list(fryDir$V1)
my.samplenames <- c("cond1r1", "cond1r2", "cond2r1", "cond2r2")
sce <- list()
sobj <- list()

for (i in 1:length(fryDir)) {
  # Load data using loadFry()
  sce[[i]] <- loadFry(fryDir[[i]], outputFormat = "scRNA", nonzero=TRUE) 
  sce[[i]] <- scater::logNormCounts(sce[[i]]) 
  # Convert to Seurat object
  sobj[[i]] <- CreateSeuratObject(counts = assay(sce[[i]], "counts"), data = assay(sce[[i]], "logcounts"))
  sobj[[i]]@meta.data$orig.ident <- my.samplenames[i]
  sobj[[i]][['RNA']] <- sobj[[i]][['RNA']]
  DefaultAssay(object = sobj[[i]]) <- "RNA"
  sobj[[i]]@project.name <- my.samplenames[i]
}
names(sobj) <- my.samplenames
```
Loading data mapped with `cellranger` directly into a `Seurat` object.
```{r}
my.data <- Read10X(data.dir = "/path/to/cellranger/output/filtered_gene_bc_matrices/")
sobj <- CreateSeuratObject(counts = my.data, project = "projectName", min.cells = 3, min.features = 200)
```
*Good saving point*: save your raw loaded data so you can go back to it and chenge filtering if needed.


## Quality control and filtering
In the quality control step we want to get rid of low quality nuclei, empty droplets with ambient RNA, and doublet nuclei in one doroplet. We start by putting together all samples for handling ease, and checking the three basic parameters we will be using for quality control: number of genes expressed per cell (`nFeature_RNA`), number of RNA molecules -- or UMIs --  per cel (`nCount_RNA`), and percentage of mitochondrial reads. The latter we will need to add manually by identifying the prefix for them in the gene names, in this case "`MT-`".
```{r}
# put the datasets together, if there are many.
if (length(sobj) > 1) {
  my.multi = T
  my.se = sobj[[1]]
  my.se = merge(my.se, sobj[2:length(sobj)], add.cell.ids = my.samplenames)
} else {my.se = sobj[[1]]}
# add mitochondrial read information
my.se[["percent.mt"]] <- PercentageFeatureSet(my.se, pattern = "^MT-")
# first look at QCs: use to set first filtering parameters
VlnPlot(my.se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", alpha = 0.1)
```
The first visualisation of the QC distributions can be used to identify a first set of threshold. Before removing them completely, we will mark them as "keep" or "not keep", so that we can plot all cells, discriminate between the ones below or above thresholds, and then refine the thresholds based on the cell-QC distributions. A good starting point for threshod are as follows:
- **Genes per cell**: a lower end of ~500, and about 10x as much for the upper limit;
- **UMIs per cell**: a lower end of ~500, or at least one UMI per gene, and about 10x as much for the upper limit
- **Mitochondrial reads**: customary agreement is a cutoff of 5%, but based on data type this could be as high as 15% (e.g. many *Drosophila* datasets)

**Important note:** quantification with `alevin-fry` will return *all* 10X droplets, including all the empty ones. For this reason, the number of detected "cells" will be much higher pre-filtering than for `cellranger`.
```{r}
keep <- rownames(subset(
	my.se, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 3)[[]])
my.se[[]]$keep <- ifelse(rownames(my.se[[]]) %in% keep, "yes", "no")
```
A violin plot of the three values we filtered on can give us a first idea of whether the thresholds we chose are good or could be immediately tweaked, based on the distribution across samples.
```{r}
p1 <- ggplot(my.se[[]], aes(x = orig.ident, y = nFeature_RNA)) +
  geom_jitter(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep)), alpha=0.1) +
  geom_jitter(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep)), alpha=0.1) +
  geom_violin(alpha=0.3) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + theme_classic()
p2 <- ggplot(my.se[[]], aes(x = orig.ident, y = nCount_RNA)) +
  geom_jitter(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep)), alpha=0.1) +
  geom_jitter(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep)), alpha=0.1) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  geom_violin(alpha=0.3) +
  #ylim(0, 35000) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + theme_classic()
p3 <- ggplot(my.se[[]], aes(x = orig.ident, y = percent.mt)) +
  geom_jitter(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep)), alpha=0.1) +
  geom_jitter(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep)), alpha=0.1) +
  geom_violin(alpha=0.3) +
  ylim(0, 25) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + theme_classic()
plot_grid(p1, p2, p3, align="h", ncol=3)
```
The first real filtering plot is a scatter plot of the number of transcript per cell to the number of genes per cell: this relationship should be somewhat linear. Low values for both indicate empty droplets, whereas high values usuall drop off the linear relationship and indicate doublets. Based on this plot, we can tweak our cut-offs.
```{r}
ggplot(my.se[[]], aes(x = nFeature_RNA, y = nCount_RNA)) +
  geom_point(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep))) +
  geom_point(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep))) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  facet_wrap(~orig.ident) +
  theme_classic()
```
Another useful scatter plot is the number of genes to percentage of mitochondrial reads. Here you expect to see an exponential decay relationship between the two values, with high percentages at low gene values, rapidly decreasing.
```{r}
ggplot(my.se[[]], aes(x = nFeature_RNA, y = percent.mt)) +
  geom_point(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep))) +
  geom_point(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep))) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  facet_wrap(~orig.ident) +
  theme_classic()
```
Once we have settles on reasonable filtering cut-off, we can effectively subset the `Seurat` object, and split it back into the individual samples, in a list.
```{r}
my.se <- subset(my.se, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 3)
my.samples <- SplitObject(my.se, split.by="orig.ident")
```
The final filtering step is to **remove doublets** that are still present after our first filtering step. Here we use `scDblFinder`, which requires the data to be a `SingleCellExperiment` object. 
```{r}
set.seed(2759)
dbl_list <- list()
for (i in seq_along(my.samples)) {
  sample <- my.samples[[i]]
  # Convert to SingleCellExperiment object
  sample_sce <- as.SingleCellExperiment(sample)
  # Calculate DoubletScore
  sample_sce$DoubletScore <- computeDoubletDensity(sample_sce)
  # Perform doublet thresholding
  sample_sce$Doublets <- doubletThresholding(data.frame(score = sample_sce$DoubletScore),
                                             method = "griffiths", returnType = "call")
  dbl_list[[i]] <- sample_sce
}
dbl_score <- function(sce) {
  plotColData(sce, x="ident", y="nCount_RNA", colour_by="DoubletScore") +
    theme_classic() + 
    theme(legend.position = "none")
}
```
Finally, we assing a "doublet" or "singlet" status to each cell, and use that to filter the `SingleCellExperiment` object, before transforming it back to a `Seurat` object.
```{r}
# filter out doublets and convert back to seurat object
sub_sce <- list()
my.sobj <- list()
for (i in seq_along(dbl_list)) {
 sub_sce <- dbl_list[[i]][, !dbl_list[[i]]$Doublets=="doublet"]
 my.sobj[[i]] <- as.Seurat(sub_sce)
 DefaultAssay(object = my.sobj[[i]]) <- "RNA"
}
names(my.sobj) <- my.samplenames
```
*Good saving point*: Save your unnormalised but filtered data for eventual changes in integration and normalisation strategies.


## Normalisation and integrating multiple datasets
There are many reasons why cells within and between samples might cluster based on gene expression, several of which are not usually interesting in the light of our biological questions. These factors include cell cycle status, batch effects, overall RNA amount, etc. So before we do any sort of integration and dimensionality reduction, we regress out these sources of uninteresting variation. What we choose to regress out depends on what we are interested in studying.
```{r}
set.seed(374)
my.cc <- cc.genes.updated.2019
my.samples <- lapply(my.samples, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
  x <- CellCycleScoring(x, s.features = my.cc$s.genes, g2m.features = my.cc$g2m.genes)
  return(x)
})
my.samples <- lapply(X = my.samples, FUN = SCTransform,
                    vars.to.regress = c("nCount_RNA","S.Score", "G2M.Score"),
                    verbose = FALSE, return.only.var.genes = F)
```
When we have multiple samples -- either from replicates of from different conditions -- we want to integrate them in order to compare their expression profiles. The integration is based on a subset of the expressed and variable genes, selected using `SelectIntegrationFeatures()`. These can be a specific amount of genes, or a proportion of the overall variable features; in this example, we take 3000 of the variable genes. These will be used as *anchors* between the different samples.
```{r}
my.features <- SelectIntegrationFeatures(object.list = my.samples, nfeatures = 3000)
my.samples <- PrepSCTIntegration(object.list = my.samples, anchor.features = my.features, verbose = FALSE)
my.anchors <- FindIntegrationAnchors(my.samples, normalization.method = "SCT", anchor.features = my.features,
                                    verbose = FALSE, dims = 1:50, k.filter = 100)
my.se <- IntegrateData(anchorset = my.anchors, normalization.method = "SCT", verbose = FALSE, preserve.order=TRUE)
```
*Good saving point*: save your integrated data before clustering, so you can play around with different dimensions and resolution.


## Dimensionality reduction and clustering
Transcriptomic datasets are *high dimensional datasets*, but our brains can only really comprehend and work in 3 of those dimensions. These steps will reduce the dimensions to 2 that we can plot, by collapsing all others baased on how much they drive variability. If this makes no sense whatsoever to you, I recommend the [PCA](https://www.youtube.com/watch?v=FgakZw6K1QQ), [UMAP](https://www.youtube.com/watch?v=eN0wFzBA4Sc), and [tSNE](https://www.youtube.com/watch?v=NEaUSP4YerM) videos from **Joshua Stramer** on **StatQuest**.

Now that we have our different samples in a single dataset, we need to find the new *highly variable features* (HVFs). These will be the genes we use for our *principle component analysis* (PCA), the first step in **dimensionality reduction**. 
```{r}
# get variable features
my.se <- FindVariableFeatures(my.se, method="sct", assay = "integrated")
# Only the genes with variability > median
my.HVF <- HVFInfo(my.se, method="sct")
my.HVF <- rownames(my.HVF)[which(my.HVF[,3] > (median(my.HVF[,3])+mad(my.HVF[,3])))]

# Run the actual PCA
my.se <- RunPCA(my.se, assay="integrated", group.by = "orig.ident", verbose=F, features=my.HVF)
# plot
DimPlot(my.se, reduction = "pca", group.by = "orig.ident")
```
Based on the PCA separation of the samples in the integrated dataset, we can assess how well the integration worked: we want the different colours to be as overlapping as possible.

PCA is great for fast linear collapsing of high dimensional data onto two axes, but it is not great for exploratory visualisation usually required in genomics. For this non-linear, stochastic methods like **t-SNE** (*t-Distributed Stochastic Neighbor Embedding*) and **UMAP** (*Uniform Manifold Approximation and Projection*) are preferrable. In this walkthrough we use UMAP, but t-SNE can be used just as well -- as long as we keep in mind the differences in their interpretation (or over-interpretation). The main and simplified difference between the two approaches can be summarised as such: t-SNE better preserves local relationships but distorts global structure, while UMAP is better at global structure fidelity at the expense of local relationships.

When running our UMAP (or t-SNE), we will need to choose the number of PC dimensions to work with, so that they are as informative as possible without becoming redundant. there are many ways to choose them, but from experience, in most comparative transcriptomics datasets the first 30 dimensions or so are always informative. To see how many more still contain strong variance, I plot them as a heatmap: when I stop seeing a strong clustering in the genes therein contained, I pull the cut-off. 
```{r}
DimHeatmap(my.se, dims=24:50, cells=500, balanced=TRUE)
```
In the example above we have good clustering of high vs low expression scores well up to dimension 50, so we run the UMAP on all of these. Since UMAP is a **stochastic** approach, theroetically it will change every time we run it anew: it is important to set a random seed (`set.seed()` at the beginning or `seed.use` in the `RunUMAP()` function) that we can call again in the future to get reproducible plots. To further check how well our integration performed, we can plot the UMAP coloured by sample.
```{r}
my.se <- RunUMAP(my.se, reduction="pca", dims=1:50, seed.use=763)
# check for weird batch effects
DimPlot(my.se, reduction="umap", group.by="orig.ident")
```
The next step that will deeply affect any downstream analysis is the **cell clustering**. At this stage we try to find biologically sensible grouping of cells based on their transcriptional profiles: this required a decent level of knowledge and understanding of our samples, i.e. how many cells types are usually in this tissue, what are the main marker genes of different statuses, etc. I highly recommend playing with the **resolution** to find the sweet spot for how many clusters are identified. Higher resolutions (>1) will lead to more, smaller clusters; smaller resolutions (up to 0.3) will lead to fewer, larger clusters. As shown in the large chunk of code below, we will also then compare the marker genes in each cluster, merging the ones that appear very similar. This is all part of the fine tuning of what is -- in the end of the day -- a completely arbitrary value.
```{r}
my.se <- FindNeighbors(my.se, reduction="pca", dims=1:50)
my.se <- FindClusters(my.se, resolution=0.5, random.seed=137)

# recalculate variable features
DefaultAssay(my.se) <- "RNA"
my.se <- NormalizeData(my.se, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
my.se <- FindVariableFeatures(my.se, assay = "RNA")
my.se <- ScaleData(my.se, features = rownames(my.se), verbose = FALSE)
my.HVF <- HVFInfo(my.se, assay = "RNA")
my.HVF <- rownames(my.HVF)[which(my.HVF[,3] > (median(my.HVF[,3])))]

# Merge pairs of clusters with less than n DEGs between them
keep.check <- T
while (keep.check == T) {
  keep.check = F
# Check the tree of clusters, to see what's the relationship between them
my.se=BuildClusterTree(my.se, dims = my.dimensions, verbose = F)
plot(Tool(object = my.se, slot = 'BuildClusterTree'))
# Check only the terminal sisters
to.check = ips::terminalSisters(my.se@tools$BuildClusterTree)
  for (i in to.check) {
    # DE between the sisters
    my.DE = FindMarkers(my.se, i[1], i[2], test.use = "MAST", #latent.vars = c("CC.Difference"),
                         min.pct = 0.25, verbose = F, assay = "RNA", features = my.HVF)
    my.DE = my.DE[which(abs(my.DE$avg_log2FC)>0.5),]
    my.lDE = length(which(my.DE$p_val_adj<0.05))

    # If less than 5, merge, and repeat
    if (my.lDE < 5) {
      my.DE = rownames(my.DE)
      cat(paste0(my.lDE, " genes differentially expressed between clusters ",i[1]," and ",i[2]," merging \n",
                  "The genes are: \n",my.DE, "\n \n"))
      my.se <- SetIdent(my.se, cells = WhichCells(my.se, idents = i[2]), value = i[1])
      keep.check = T
    }
  }
}

# renumber starting from 1 (of course optional)
my.ID <- factor(Idents(my.se),levels= levels(Idents(my.se))[ base::order(as.numeric(levels(Idents(my.se))))])
levels(my.ID) <- 1:length(levels(my.ID))
my.se[["seurat_clusters"]] <- my.ID
Idents(my.se) <- "seurat_clusters"

# Plot and check clustering
DimPlot(my.se, reduction = "umap", group.by = "seurat_clusters", label=TRUE, raster=FALSE)
```
Once we are happy with our clustering, it is always good to do a QC check on each cluster, to identify some cells that might be clustering together just because they share the same low quality flags.
```{r}
# make sure you are pulling out "nFeature_RNA", "nCount_RNA", and "percent.mt"
qc.df <- melt(my.se[[]][,c(2,3,4,15)])
ggplot(qc.df, aes(x=seurat_clusters, y=value, fill=seurat_clusters)) +
	geom_point(position = position_jitter(seed = 1, width = 0.2), alpha=0.1, size=0.5) + 
	geom_violin(alpha=0.7) +
	facet_wrap(~variable, ncol=3, scales = "free") +
	theme_classic() +
	scale_x_discrete(guide = guide_axis(angle = 45))
```
*Good saving point*: Now that you have your processed and clustered sample, you can save it before annotating cell types.


## Marker genes and cell type annotation
To annotate the clusters as cell types, we need to know which genes are driving their clustering, i.e. which are the marker genes per cluster. The `Seurat` function `FindAllMarkers()` does exactly that; in this example, we are calculating marker genes from our HVFs that are expressed in at least 25% of cells in a cluster (`min.pct = 0.25`), that show a log fold-change of at least 0.5 (`logfc.threshold = 0.5`), and that have a *p*-value below 0.05 (`return.thresh = 0.05`). All these values can be played with if we want to increase or decrease stringency.

To help us inspect the markers of each cluster, I recommend pulling out the top 20 and top 3 markers per cluster: the former can be used to look up unfamiliar genes and where they might be expressed in literature, the latter can be plotted to already annotate easily distinguishible cell types.
```{r}
semarkers <- FindAllMarkers(my.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,
                           assay = "RNA", features = my.HVF, return.thresh = 0.05, verbose = F)
# embed this information in the object
my.se@misc$semarkers <- semarkers
# Take only the top for manual inspection
top20 <- semarkers %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.01) %>% 
  arrange(cluster, desc(avg_log2FC)) %>% 
  group_by(cluster) %>% top_n(20, avg_log2FC) 
top3 <- semarkers %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.01) %>% 
  arrange(cluster, desc(avg_log2FC)) %>% 
  group_by(cluster) %>% top_n(3, avg_log2FC)
```
Look at the top 3 per cluster: this will already give you an idea of the identity of many clusters based on their biology.
```{r}
DotPlot(my.se, features = unique(top3$gene), dot.scale = 6) +
    RotatedAxis() + scale_y_discrete(limits=rev) + scale_colour_viridis()
```
For others you might want to look into the expression of known markers.
```{r}
DotPlot(my.se, features = toupper(c(
	# Dermal condensate
	"Sox2", "Twist2", "Dkk1", "Sox9","Trps1",
	# Placodes
    "Edar", "Ctnnb1", "Wnt2", "Wnt7b", "Wnt10a", "Wnt10b",
	# Keratinocytes
	"Krt1", "Krt10", "Krt14",
	# Other fibroblasts
    "Pdgfra", "Fbln1"),
	dot.scale = 6) + RotatedAxis() + scale_y_discrete(limits=rev) + scale_colour_viridis()
```
*Good saving point*:
