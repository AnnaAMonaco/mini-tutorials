## Differential expression analysis using DESeq2 on RNA-seq data
Below are the packages needed for this tutorial.
```
# for DEG analysis
library(tximport)
library(tximportData)
library(tidyverse)
library(DESeq2)
library(edgeR)
# for data wrangling
library(dplyr)
library(stringr)
library(readr)
library(reshape2)
# for plotting
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library(hrbrthemes)
```
### Load count data
Start by loading the data quantified with `salmon`; you will need a tab-separated metadata file containing the minimal columns `Run`, `name`, `replicate`, and any other categorical variable that is relevant for your comparison. `Run` must be the name of the sample-specific directory that comes as output from `salmon`, and in the main directory `dir` where your metadata file also lives. The transcript-to-gene (`tx2g`) file is another tab-separated file where one column is the name of the transcript, and the other the name of the gene.
```
dir <- "/path/to/salmon_quant/"
samples <- read.table(file.path(dir,"metadata.txt"), header=TRUE)
# tell R where to find the files
files <- file.path(dir, samples$Run, "quant.sf")
names(files) <- samples$name
# load transcript to gene files
# this tells txi which transcripts correspond to which gene, and collapses counts
tx2g <- read.delim("/path/to/txp2gene.tsv", header=FALSE, sep="\t")
# import sample quantifications
txi <- tximport(files, type="salmon", tx2gene=tx2g)
```

### Make DDS object and check quality
You want to have all samples you will be comparing in one `dds` object, this is important for the normalisation step. In the `dds` generation, `design` tells you what you are comparing; if you have multiple variables, you can combine them as `~var1+var1` or check their interaction with `~var1:var2`, and the most important variable for you comparison always goes last.
```
dds <- DESeqDataSetFromTximport(
          txi,
          colData=samples,
          design=~tissue+genotype)
colnames(dds) <- samples$name

# it's standard practice to pre-filter out genes with less than 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```
Before differential expression analysis, we need to compare how similar the samples are to check replicate quality. We do this with two methods: Euclidean distance and principal component analysis (PCA).
```
# variance stabilised transformation for normalisation
vsd <- vst(dds, blind = FALSE)
# calculate Euclidean distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
ord <- hclust(dist(t(assay(vsd))))$order
# store plot
p1 <- ggplot(melt(sampleDistMatrix), aes(Var1, Var2, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  scale_x_discrete(limits=melt(sampleDistMatrix)$Var1[ord], guide = guide_axis(angle = 90)) +
  scale_y_discrete(limits=melt(sampleDistMatrix)$Var1[ord]) +
  theme_ipsum()
# calculate PCA
pcaData <- plotPCA(vsd, intgroup=c("tissue","genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# store plot
p2 <- ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic()
# generate plot
plot_grid(p1,p2,align="h")
#save plot
ggsave("img/Distance-PCAplot.png")
```
### Calculate differentially expressed genes (DEGs)
The `dds` object needs to be normalised using `DESeq2` before running contrasts. When running the contrast with `results()`, the order of the variables to intersect is important: the first one is positive LFC, the second negative LFC. Multiple contrasts can be stored in a list and processed in parallel for subsequent steps.
```
# normalisation step
dds <- DESeq(dds)
res <- list()
# alpha and lfcThreshold are the p value and LFC cutoff, can be adjusted
res[[1]] <- results(dds, contrast=c("genotype", "GT1", "GT2"), alpha=0.05, lfcThreshold=1, pAdjustMethod="BH")
res[[2]] <- results(dds, contrast=c("genotype", "GT1", "GT3"), alpha=0.05, lfcThreshold=1, pAdjustMethod="BH")

# you an check the summary of number of DEGs using summary for each contrast
for (i in seq_along(res)) {
  summary(res[[i]])
}

# this function filters DEGs for LFC of min 0.5 and p val below 0.05
filter_dataframe <- function(df) {
  df %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::arrange(desc(log2FoldChange))
}
f.res <- lapply(res, filter_dataframe)

# save the rds objects to be loaded in the future
saveRDS(res, "data/Rdata/Exp1-DEGs-lfc1-results.rds")
saveRDS(f.res, "data/Rdata/Exp1-DEGs-lfc1-results-filtered.rds")
```

Volcano plots
```
t.res <- list()
# this function assigns info we can use for aesthetic values later
add_thresh <- function(df) {
  df %>%
    data.frame() %>%
    mutate(threshold = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "high",
        padj < 0.05 & log2FoldChange < -1 ~ "low",
        padj >= 0.05 ~ "none"))
}
t.res <- lapply(res, add_thresh)
saveRDS(t.res, "data/Rdata/Aa-DEGs-lfc1-results-threshold.rds")

Mm.t.res <- list()
Mm.t.res <- lapply(Mm.res, add_thresh)
saveRDS(Mm.t.res, "data/Rdata/Mm1214-DEGs-lfc1-results-threshold.rds")

goi <- list()
for (i in seq_along(f.res)) {
  goi[[i]] <- unique(c((f.res[[i]] %>% arrange(padj) %>% head(n=10))$gene,
    (f.res[[i]] %>% arrange(log2FoldChange) %>% head(n=10))$gene,
    (f.res[[i]] %>% arrange(log2FoldChange) %>% tail(n=10))$gene))
}
# example of volcano plot
p1 <- ggplot(t.res[[1]], aes(x = log2FoldChange, y = -log10(padj), colour = thre
shold, alpha=threshold)) +
  geom_point() +
  scale_alpha_manual(name="threshold", values = c("high" = 1, "low" = 1, "none" 
= 0.1)) +
  scale_colour_manual(name="threshold", values = c("high" = "#E91E63", "low" = "
#EF9A9A", "none" = "grey78")) +
  ggtitle("Aa E22.5 Do vs Ve") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(label=ifelse(rownames(t.res[[1]]) %in% goi[[1]], rownames(t.res[[1]]), ""), size=3, max.overlaps=30, colour = "black") +
  theme_classic() +
  theme(legend.position = "none")
  ggsave("img/MmAa-E14DOvsVE-volcano.pdf")
```

### other plots
Relative expression boxplot
```
Mm.df <- counts(Mm.dds, normalized=TRUE)["Tbx15",] %>% melt %>%
  mutate(
    skin=case_when(
        grepl("Do", rownames(.)) ~ "Dorsal",
        grepl("Ve", rownames(.)) ~ "Ventral"
      ),
    time=case_when(
        grepl("Do14", rownames(.)) ~ "2_placode",
        grepl("Do12", rownames(.)) ~ "1_preplacode",
        grepl("Ve14", rownames(.)) ~ "1_preplacode"
        )) 
p <- list()
p[[1]]<- ggplot(Aa.df, aes(x=skin, y=value, fill=skin, colour=skin)) + 
  geom_boxplot(alpha=0.4) +
  scale_fill_manual(values = c("Dorsal" = "#E91E63", "Ventral" = "#EF9A9A")) +
  scale_color_manual(values = c("Dorsal" = "#E91E63", "Ventral" = "#EF9A9A")) +
  scale_y_log10(limits=c(10,100000)) +
  facet_wrap(~ time) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
p[[2]]<- ggplot(Mm.df, aes(x=skin, y=value, fill=skin, colour=skin)) + 
  geom_boxplot(alpha=0.4) +
  scale_fill_manual(values = c("Dorsal" = "#FBC02D", "Ventral" = "#FFECB3")) +
  scale_color_manual(values = c("Dorsal" = "#FBC02D", "Ventral" = "#FFECB3")) +
  scale_y_log10(limits=c(10,100000)) +
  facet_wrap(~ time) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
plot_grid(plotlist=p, align="v")
ggsave("img/Tbx15-counts.pdf")
```
Differential expression heatmaps
```
```
