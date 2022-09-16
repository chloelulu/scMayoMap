<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMayoMap

<!-- badges: start -->
<!-- badges: end -->

The goal of scMayoMap is to annotate cell clusters for single cell data.

# 1. Installation

``` {r, message=FALSE, warnings = F, results=FALSE}
devtools::install_github("chloelulu/scMayoMap")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

Library required packages
```{r, message=FALSE, warning=FALSE}
library(scMayoMap)
```

After you load the package, the demo “data” and integrated “scMayoMapDatabase”
will be automatically loaded

``` {r, message=FALSE, warning=FALSE}
head(data)
```

scMayoMapDatabase is the database integrated in the package. User can also supply their own database, should be in the same format as scMayoMapDatabase

``` {r, message=FALSE, warning=FALSE}
head(scMayoMapDatabase[,c(1:4)])
```


# 2. Run scMayoMap

To annotate a scRNA-seq data analysis output generated by [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers)
You can specify the tissue type (i.e.,tissue = 'muscle') or not (tissue = NULL, default). But we do suggest to specify the tissue type, which will be more accurate. 

```{r, message=FALSE, warning=FALSE}
obj <- scMayoMap(data = data, tissue = 'muscle')
```

User can also define their own database as input. For example, if you want to use your own marker pool, for [example](https://github.com/chloelulu/scMayoMap/blob/main/data/demo.marker.Rdata).
The database should be in the same format as scMayoMapDatabase. Presence of a gene in a celltype will be 1, otherwise 0. column names should contain 'gene' column. Currently, our database only supports gene symbols. If you have a data with Ensembl ID, please refer [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) or other methods to convert the gene names. 
```{r, message=FALSE, warning=FALSE}
load('data/demo.marker.Rdata')
db <- tidyr::spread(demo.marker, key = c('celltype'), value = 'value')
db[is.na(db)] <- 0
obj <- scMayoMap(data = data, database = database)
```


You can extract the result summary table and export it according to your preference.
```{r, message=FALSE, warning=FALSE}
res <- obj$res
head(res)
```

You can extract the marker genes table used for prediction. See descriptions of each column in scMayoMap R Documentation.
```{r, message=FALSE, warning=FALSE}
markers <- obj$markers
head(markers)
```

scMayoMap also supplies a plot function to help you visualize the result. You can supply your own directory(i.e., directory = '~/Desktop/'), otherwise plot will be automatically saved to the current working directory (default is "directory = NULL").
```{r, fig.retina = 4, fig.width= 7, fig.height=5, results=FALSE, message=FALSE, warning=FALSE}
plt <- scMayoMap.plot(scMayoMap.object = obj, directory = '~/Desktop/', width = 8, height = 6)
```

# 3. Using scMayoMap start from Seurat 

Below code illustrates how to use scMayoMap if you start with a Seurat analysis. We will use a demo dataset from [SeuratData](https://github.com/satijalab/seurat-data). You need to install package [Seurat](https://satijalab.org/seurat/articles/install.html) and [SeuratData](https://github.com/satijalab/seurat-data) before try below example. However, if you have a large dataset, you may refer to [multi-thread](https://github.com/satijalab/seurat/issues/1865) to speed up [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) procedure.
```{r,fig.retina = 4, fig.width= 7, fig.height=5, results=FALSE, message=FALSE, warning=FALSE}
library(SeuratData)
library(Seurat)
seurat.obj <- LoadData(ds = "pbmc3k") 
seurat.obj <- NormalizeData(object = seurat.obj)
seurat.obj <- FindVariableFeatures(object = seurat.obj)
seurat.obj <- ScaleData(object = seurat.obj)
seurat.obj <- RunPCA(object = seurat.obj)
seurat.obj <- FindNeighbors(object = seurat.obj)
seurat.obj <- FindClusters(object = seurat.obj)
seurat.markers <- FindAllMarkers(seurat.obj, method = 'MAST')
scMayoMap.obj <- scMayoMap(data = seurat.markers, database=scMayoMapDatabase, tissue = 'blood')
plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
```

You can also extract the markers mapped by scMayoMap and making a dotplot. 
Below code visualize the genes mapped to cluster 0 (predicted cell type Granulocyte)
```{r,fig.retina = 4, fig.width= 8, fig.height=5, results=FALSE, message=FALSE, warning=FALSE}
gns <- scMayoMap.obj$markers$genes[scMayoMap.obj$markers$cluster==0 & scMayoMap.obj$markers$celltype=='Granulocyte']
gns <- strsplit(gns, ',')[[1]]
DotPlot(seurat.obj, features = gns)
```

# 4. References

  - Lu Yang, Yan Er Ng, Haipeng Sun, Ying Li, Nathan LeBrasseur, Jun Chen, Xu Zhang. Single-cell Mayo Map (scMayoMap), an easy-to-use tool for mapping cell types in single-cell RNA-sequencing data analysis. Submitted.

# scMayoMap
