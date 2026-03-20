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
Loading data mapped with `cellranger` directly into a `Surat` object.
```{r}
```

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
The first visualisation of the QC distributions can be used to identify a first set of threshold, Before removing them completely, we will mark them as "keep" or "not keep", so thata we acn plot all cells, discriminate brtween the ones below or above thresholds, and then refine the thresholds based on the cell-QC distributions. A good starting point for threshod are as follows:
- Genes per cell: a lower end of ~500, and about 10x as much for the upper limit;
- UMIs per cell: a lower end of ~500, or at least one UMI per gene, and about 10x as much for the upper limit
- Mitochondrial reads: customary agreement is a cutoff of 5%, but based on data type this could be as high as 15% (e.g. many *Drosophila* datasets)

**Important note:** quantification with `alevin-fry` will return *all* 10X droplets, including all the empty ones. For this reason, the number of detected "cells" will be much higher pre-filtering than for `cellranger`.
```{r}
keep <- rownames(subset(
	my.se, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 3)[[]])
my.se[[]]$keep <- ifelse(rownames(my.se[[]]) %in% keep, "yes", "no")
```
```{r}
ggplot(my.se[[]], aes(x = nFeature_RNA, y = nCount_RNA)) +
  geom_point(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep))) +
  geom_point(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep))) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  facet_wrap(~orig.ident) +
  theme_classic()
```

```{r}
ggplot(my.se[[]], aes(x = nFeature_RNA, y = percent.mt)) +
  geom_point(data = subset(my.se[[]], keep == "yes"), aes(col = as.factor(keep))) +
  geom_point(data = subset(my.se[[]], keep == "no"), aes(col = as.factor(keep))) +
  scale_color_manual(name = "keep", values = c("yes" = "grey18", "no" = "brown1")) +
  facet_wrap(~orig.ident) +
  theme_classic()
```

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

```{r}
my.se <- subset(my.se, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 3)
my.samples <- SplitObject(my.se, split.by="orig.ident")
```
Remove doublets
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
plots <- lapply(dbl_list, dbl_score)
plot_grid(plotlist=plots, align="h", ncol=2)
```
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


## Normalisation and integrating multiple datasets
We regress out sources of uninteresting variation, like cell cycle status and overall RNA amount 
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
We then integrate out datasets based on a subset of the expresse genes, selected using `SelectIntegrationFeatures()`. These can be a specific amount of genes, or a proportion of the overall variable features.
```{r}
my.features <- SelectIntegrationFeatures(object.list = my.samples, nfeatures = 3000)
options(future.globals.maxSize = 30 * 1024^3)  # 30 GiB
my.samples <- PrepSCTIntegration(object.list = my.samples, anchor.features = my.features, verbose = FALSE)
my.anchors <- FindIntegrationAnchors(my.samples, normalization.method = "SCT", anchor.features = my.features,
                                    verbose = FALSE, dims = 1:50, k.filter = 100)
my.se <- IntegrateData(anchorset = my.anchors, normalization.method = "SCT", verbose = FALSE, preserve.order=TRUE)
```

## Dimensionality reduction and clustering
Find variable features and run PCA
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
Based on PCA separation of samples you can assess how well the integration worked.

Check dimension heatmap to identify cutoff of dimensions to use.
```{r}
DimHeatmap(my.se, dims=24:50, cells=500, balanced=TRUE)
```
Dimensionality reduction with UMAP
```{r}
my.se <- RunUMAP(my.se, reduction = "pca", dims = 1:50, seed.use=763)
# check for weird batch effects
DimPlot(my.se, reduction = "umap", group.by = "orig.ident")

```
Clustering: play around with resolution to find sweet spot based on biological knowldge or downstream goals.
```{r}
my.se <- FindNeighbors(my.se, reduction = "pca", dims = 1:50)
my.se <- FindClusters(my.se, resolution = 0.5, random.seed = 137)

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
Quick QC check per cluster
```{r}
qc.df <- melt(my.se[[]][,c(2,3,4,15)])
ggplot(qc.df, aes(x=seurat_clusters, y=value, fill=seurat_clusters)) +
	geom_point(position = position_jitter(seed = 1, width = 0.2), alpha=0.1, size=0.5) + 
	geom_violin(alpha=0.7) +
	facet_wrap(~variable, ncol=3, scales = "free") +
	theme_classic() +
	scale_x_discrete(guide = guide_axis(angle = 45))
```

## Marker genes and cell type annotation
Start by finding all the marker genes per cluster.
```{r}
semarkers <- FindAllMarkers(my.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold =  0.5,
                           assay = "RNA", features = my.HVF, return.thresh = 0.05, verbose = F)
my.se@misc$semarkers <- semarkers
# Take only the top
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
