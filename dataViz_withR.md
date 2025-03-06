# Basic visualisation in R

## Requirements for this mini-tutorial


## Basic R vs `ggplot2` 

### Aesthetics in `ggplot2`

#### Working with colours

why palettes are important and how to make sure they have high contrast and are colourblind friendly.
some nice pre-existing palettes (viridis, cartocolors, wesanderson)
fill vs colour 
gradient vs factorial
consistent colour code

##### AVOID: Bidirectional colour scales for unidirectional data

#### Size (and shape) sometimes matters

how size can add an additional scalar measurement on top of colour
how shape can add another categorical variable on top of colour


## Main types of plots you are likely to use

### Dot plot

#### Scatter plot

##### Lion's Q: How to best visualise bulk RNA-seq?
My opinion: MA plots, volcano plots, expression scatterplots.

#### Bubble plot

#### Linking point with a line

### Bar plot

#### Proportion stacked bar plots instead of pie charts
Humans do not intuitively understand angle-based differences

##### AVOID: Bar plots for means separation

##### AVOID: Confusing position-based viz with length-based viz

##### AVOID: Not checking the data range at each factor level

### Box plot

### Violin plot

##### AVOID: Violin plots for small sample sizes

### Heatmap

##### AVOID: Unclustered heatmaps

##### AVOID: Heatmaps without checking for outliers


## Composite plots

### Faceting avoids overwhelming "meadows"

### Violin plot with boxplot inside
You get information on both mean and distribution

### Adding data points to other plots

### Fitting smoothened lines to plots

### Adding statistical information
Error bars, mean, median, etc.







