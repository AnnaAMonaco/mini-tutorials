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

## Quality control and filtering

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

```{r}
VlnPlot(my.se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", alpha = 0.1)
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


## Integrating multiple datasets
## Normalisation
## Dimensionality reduction and clustering
## Marker genes and cell type annotation
