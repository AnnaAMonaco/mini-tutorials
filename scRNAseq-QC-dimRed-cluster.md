Here you will learn the basics of analysing scRNA-seq datasets in R. I will cover the following steps:
- Loading mapped data from salmon alevin and from cellranger
- Quality control and filtering
- Integrating multiple datasets
- Normalisation
- Dimensionality reduction and clustering
- Marker genes and cell type annotation

These are the libraries you will need to have installed and loaded:
```{r}
library(fishpond)
library(tximport)
library(tximeta)
library(SingleCellExperiment)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(scran)
library(scater)
library(scDblFinder)
library(BiocSingular)
library(scales)
```

## Load data
If you quantified with `alevin-fry` you can load the dataset with `loadFry()`. The resulting object is a `SingleCellExperiment` object, so it will also need to be transformed into a `SeuratObject` for downstream processing.
```{r}
fryDir <- read.delim("samples.txt", header=FALSE)
fryDir <- as.list(fryDir$V1)
my.samplenames <- c("AaDoPcR1", "AaDoPcR2", "AaDoPcR3", "AaVePcR1", "AaVePcR2")
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
## Integrating multiple datasets
## Normalisation
## Dimensionality reduction and clustering
## Marker genes and cell type annotation
