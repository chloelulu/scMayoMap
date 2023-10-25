<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMayoMap

<!-- badges: start -->
<!-- badges: end -->

scMayoMap aims to annotate cell clusters in single-cell data, facilitating the analysis and interpretation of single-cell RNA-sequencing (scRNA-seq) data.

# 1. Installation
Install the development version of scMapping using the following commands:

``` {r, message=FALSE, warnings = F, results=FALSE}
devtools::install_github("chloelulu/scMayoMap")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

Load the necessary packages:
```{r, message=FALSE, warning=FALSE}
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
sapply(pkgs, require, character.only = TRUE)
library(scMayoMap)
```

Upon loading scMayoMap, the demo "data" and integrated "scMayoMapDatabase" will automatically be available.

``` {r, message=FALSE, warning=FALSE}
head(data)
```

scMayoMapDatabase is included within the package. Users can also utilize their marker pool, ensuring they match the format of scMappingDatabase.

``` {r, message=FALSE, warning=FALSE}
head(scMayoMapDatabase[,c(1:4)])
```


# 2. Run scMayoMap

To annotate a scRNA-seq data analysis output generated by [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers)
You can specify the tissue type (i.e.,tissue = 'muscle') or not (tissue = NULL, default). But we do suggest to specify the tissue type, which will be more accurate. Additionally, please make sure your data is in the same format as demo data. 

```{r, message=FALSE, warning=FALSE}
obj <- scMayoMap(data = data, tissue = 'muscle')
```

User can also define their own database as input. For [example](https://github.com/chloelulu/scMayoMap/blob/main/data/demodata.rda), if you want to use your own marker pool. The database should be in the same format as scMayoMapDatabase. Presence of a gene in a celltype will be 1, otherwise 0. column names should contain 'gene' column. Currently, our database only supports gene symbols. If you have a data with Ensembl ID, please refer [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) or other methods to convert the gene names. 
```{r, message=FALSE, warning=FALSE}
head(demodata)
db <- tidyr::spread(demodata, key = c('celltype'), value = 'value')
db[is.na(db)] <- 0
obj <- scMayoMap(data = data, database = db)
```


You can retrieve the result summary table and export it in a format of your choosing.
```{r, message=FALSE, warning=FALSE}
res <- obj$res
head(res)
```

You can retrieve the table of marker genes utilized for predictions. Refer to the scMayoMap R Documentation for detailed descriptions of each column.
```{r, message=FALSE, warning=FALSE}
markers <- obj$markers
head(markers)
```

scMayoMap provides a plotting function to aid in visualizing the results. You can specify a directory where the plot will be saved (e.g., directory = '~/Desktop/'); if none is provided, the plot will be automatically saved to the current working directory (default is "directory = NULL"). Users should adjust the "width" and "height" parameters to suit their preferred figure size.
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

You can also retrieve the markers identified by scMayoMap to create a dot plot. The following code illustrates how to visualize the genes mapped to cluster 0, where the predicted cell type is Granulocyte.
```{r,fig.retina = 4, fig.width= 8, fig.height=5, results=FALSE, message=FALSE, warning=FALSE}
gns <- scMayoMap.obj$markers$genes[scMayoMap.obj$markers$cluster==0 & scMayoMap.obj$markers$celltype=='Granulocyte']
gns <- strsplit(gns, ',')[[1]]
DotPlot(seurat.obj, features = gns)
```

# 4. References
If you use scMayoMap, please cite: 
  - Lu Yang, Yan Er Ng, Haipeng Sun, Ying Li, Nathan LeBrasseur, Jun Chen, Xu Zhang. Single-cell Mayo Map (scMayoMap), an easy-to-use tool for mapping cell types in single-cell RNA-sequencing data analysis. Submitted.

