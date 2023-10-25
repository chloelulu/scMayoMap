
#' @title Calculate the cumulative variance of a vector
#' @param x a vector
#' @return {v}{a vector of cumulative variance}
#' @rdname cumvar
#' @import stats
cumvar <- function(x){
  x <- x - x[sample.int(length(x), 1)]
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  return(v)
}


utils::globalVariables(c("celltype", "cluster", "cell","score", "scMayoMapDatabase"))

#' An illustrative user-supplied marker pool
#' This dataset serves as a demonstration and includes a collection of user-specified marker pools.
#' @docType data
#' @usage data(demodata)
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{gene}{Names of marker genes}
#'   \item{celltype}{Names of cell types}
#'   \item{value}{A column populated with the value of 1, signifying the presence of the corresponding marker gene for the respective cell type}
#' }
#' @examples
#' # If 'demodata' is your marker pool, use the code below to modify the marker
#' # pool into the database format required by scMayoMap:
#' library(scMayoMap)
#' data(demodata)
#' db <- tidyr::spread(demodata, key = c('celltype'), value = 'value')
#' db[is.na(db)] <- 0
#' obj <- scMayoMap(data = data, database = db)
"demodata"


#' scMayoMapDatabase
#' The default databse for scMayoMap. A data.frame lists all marker genes for each cell type in a variety of tissues.
#' @docType data
#' @usage data(scMayoMapDatabase)
#' @examples
#' data(scMayoMapDatabase)
"scMayoMapDatabase"


#' @title scMayoMap
#' @description This function is crafted to annotate cell clusters in single-cell RNA-Seq data and is particularly tailored to enhance the downstream analysis of the Seurat workflow. Additionally, it is versatile enough to accommodate the downstream analysis of other tools, provided the input format is compatible.
#' @param data A data frame where rows represent genes and columns include "avg_log2FC", "pct.1", "pct.2", and "gene". The "avg_log2FC" column indicates the log fold-change of average expression between each cell cluster and all others. Columns "pct.1" and "pct.2" display the percentage of cells in which the feature is detected in each cell cluster and all others, respectively. The "gene" column contains the gene names. For more details, refer to outputs of \link[Seurat]{FindAllMarkers}.
#' @param database data frame. If NULL (default), the function will automatically use the default scMayoMapDatabase. Users can also provide their own set of markers, ensuring that they are formatted similarly to scMayoMapDatabase, see \link[scMayoMap]{scMayoMapDatabase} and \link[scMayoMap]{demodata}.
#' @param tissue The tissue name to be specified.  If NULL, no tissue will be specified. Currently, scMayoMapDatabase includes the following tissues:
#'                "adipose tissue"; "bladder"; "blood"; "bone"; "bone marrow";
#'                "brain";"breast"; "embryo"; "eye"; "gastrointestinal tract";
#'                "heart";"kidney";"liver";"lung";"mammary gland";"muscle";
#'                "other";"ovary";"pancreas";"placenta";"prostate";
#'                "skin";"spleen";"stomach";"testis";"thymus";"tooth";"uterus".
#' @param padj.cutoff A numerical value between 0 and 1 indicating the significance level for keep the differential genes as input for scMayoMap. Default is 0.05.
#' @param pct.cutoff A numerical value between 0 and 1. Genes with pct.1 greater than this number will be remained for later analysis.
#' @return A list with the elements
#' \item{res}{A data frame that includes the cluster ("cluster") associated with the top N potential cell types ("celltype"), along with the estimated scores "ES.raw", standard errors "ES.se", normalized estimated scores "ES.norm", and markers matched for the predicted cell type.}
#' \item{markers}{A data frame containing the scores of cell types for each cluster. Row: cluster, Column: annotated cell type.}
#' \item{tissue}{The tissue names specified in the analysis.}
#' \item{annotation.norm}{A data frame containing the scores of cell types for each cluster. Row: cluster, Column: annotated cell type.}
#' \item{annotation.mean}{A data frame containing the mean scores of cell types for each cluster. Row: cluster, Column: annotated cell type.}
#' @rdname scMayoMap
#' @importFrom dplyr summarise
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom tibble column_to_rownames
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @export scMayoMap
#' @author Lu Yang and Xu Zhang
#' @references Lu Yang, Yan Er Ng, Haipeng Sun, Ying Li, Lucas C. S. Chini, Nathan K. LeBrasseur, Jun Chen & Xu Zhang. Single-cell Mayo Map (scMayoMap): an easy-to-use tool for cell type annotation in single-cell RNA-sequencing data analysis. BMC Biol 21, 223 (2023). https://doi.org/10.1186/s12915-023-01728-6.
#' @examples
#' pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
#' sapply(pkgs, require, character.only = TRUE)
#' library(scMayoMap)
#' data(data)
#' obj <- scMayoMap(data = data, tissue = 'muscle')

scMayoMap <- function(data, database=NULL, padj.cutoff = 0.05, pct.cutoff = 0.25, tissue = NULL){

  if(is.null(database)){
    database = scMayoMapDatabase
  }else{
    cat('User defined marker pool will be used!\n')
  }

  ## Define tissue or not
  if(is.null(tissue)){
    db <- database[,!(colnames(database) %in% 'tissue')]
  }else{
    db <- (database[database$tissue==tissue,,drop =F])[,!(colnames(database) %in% 'tissue')]
  }
  data1 <- data
  colnames(data1)[grep('log',colnames(data1))] = "avg_log2FC"

  ## Filter gene by cutoffs
  data1 <- data1[data1$p_val_adj <= padj.cutoff & data1$pct.1 >= pct.cutoff,]
  data1$gene <- toupper(data1$gene)

  ## Match the database
  min.pct2 <- min(data1$pct.2[data1$pct.2>0])
  data1$Score <- ifelse(data1$pct.2 !=0, (2^data1$avg_log2FC * data1$pct.1)/data1$pct.2,(2^data1$avg_log2FC * data1$pct.1)/min.pct2)
  merged <- na.omit(merge(x=data1[,c('cluster','gene','Score')], y=db, all.y = TRUE))
  merged[,!(colnames(merged) %in% c('gene','Score','cluster'))] <- sapply(merged[,!(colnames(merged) %in% c('gene','Score','cluster'))], function(x) merged$Score * x )
  merged1 <- merged[,!(colnames(merged) %in% c('gene','Score')), drop =F]

  ## Calculate the mean value
  annotation.mean <- as.data.frame(summarise(group_by(merged1, cluster),across(everything(), mean)))
  annotation.se <- as.data.frame(summarise(group_by(merged1, cluster), across(everything(), function(x) sd(x)/sqrt(length(x)))))
  rownames(annotation.mean) <- annotation.mean$cluster
  rownames(annotation.se) <- annotation.se$cluster
  annotation.mean <- annotation.mean[,!(colnames(annotation.mean) %in% 'cluster')]
  annotation.se <- annotation.se[,!(colnames(annotation.se) %in% 'cluster')]

  ## Renormalize
  annotation.norm <- annotation.mean / rowSums(annotation.mean)
  if(sum(is.na(annotation.norm)) > 0){cat('Caution: some clusters can not match db!')}
  annotation.norm[is.na(annotation.norm)] <- 0

  ## Summarize the top.n possible cell types
  top.n <- sapply(1:nrow(annotation.norm),  function (i) {
    x <- unlist(annotation.norm[i, ])
    y <- unlist(annotation.se[i, ])
    z <- unlist(annotation.mean[i, ])
    x <- sort(x,decreasing = T)
    xx <- cumvar(x)
    xx[1] <- 0
    n <- unname(which.max(diff(xx)))
    nms <- paste0(names(rev(x[1:n])), collapse = '; ')

    scores.norm <- paste0(round(rev(x[1:n]),2), collapse = '; ')
    y <- y[names(rev(x[1:n]))]
    SEs <- paste0(round(y[names(rev(x[1:n]))],2), collapse = '; ')
    scores.raw <- paste0(round(z[names(rev(x[1:n]))],2), collapse = '; ')
    return(c(nms, scores.norm, scores.raw, SEs))
  })

  top.n <- cbind.data.frame(t(top.n))
  top.n$cluster <- rownames(annotation.norm)
  colnames(top.n) <- c('celltype','ES.norm','ES.raw','ES.se','cluster')
  rownames(top.n) <- NULL

  ## save marker genes for mapped cell type for Dotplot
  mk.lst <- list()
  for(k in top.n$cluster){
    clt <- top.n$celltype[top.n$cluster==k]
    clt <- strsplit(clt,'; ')[[1]]
    scores <- strsplit(top.n$ES.norm[top.n$cluster==k],'; ')[[1]]
    tmp <- c()
    for(clt.sub in 1:length(clt)){
      xx <- merged[merged$cluster==k,][,c('gene',clt[clt.sub])]
      xx <- xx[order(-xx[,2]),]
      genes <- xx$gene[which(xx[,2]>0)]
      sc <- scores[clt.sub]
      gg <- paste0(genes,collapse = ',')
      sc.gg <- c(clt[clt.sub], sc, gg)
      tmp <- rbind(tmp,sc.gg)
    }
    tmp <- as.data.frame(tmp)
    tmp$cluster <- k
    colnames(tmp) <- c('celltype','score','genes','cluster')
    mk.lst[[k]] <- tmp
  }
  markers <- do.call(rbind.data.frame, mk.lst)
  rownames(markers) <- NULL
  markers <- markers[,c('cluster','celltype','score','genes')]

  if(!is.null(tissue)){
    top.n$celltype <- sapply(1:length(top.n$celltype), function(i){
      tmp <- paste0(gsub('.*:','',strsplit(unlist(top.n$celltype[i]),'; ')[[1]]),collapse = '; ')
    })
    markers$celltype <- gsub('.*:','',markers$celltype)
  }

  return(list(res = top.n, markers = markers, tissue = tissue, annotation.norm = annotation.norm, annotation.mean = annotation.mean))
}




#' @title Visualize scMayoMap Annotation Results
#' @description This function creates a dot plot visualizing the annotated cell types for cell clusters, utilizing the output from scMayoMap.
#' @param scMayoMap.object The result object returned from the function \code{scMayoMap}.
#' @param directory A character string indicating the directory where the figures will be saved. By default (NULL), it is set to the current working directory.
#' @param width Specifies the width of the graphics region in inches. Refer to the R function \code{pdf}.
#' @param height Specifies the height of the graphics region in inches. Refer to the R function \code{pdf}.
#' @return
#' \item{annotation.plot}{A ggplot object}
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom tidyr unnest
#' @importFrom tidyr spread
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom utils globalVariables
#' @rdname scMayoMap.plot
#' @export scMayoMap.plot
#' @author Lu Yang
#' @references Lu Yang, Yan Er Ng, Haipeng Sun, Ying Li, Lucas C. S. Chini, Nathan K. LeBrasseur, Jun Chen & Xu Zhang. Single-cell Mayo Map (scMayoMap): an easy-to-use tool for cell type annotation in single-cell RNA-sequencing data analysis. BMC Biol 21, 223 (2023). https://doi.org/10.1186/s12915-023-01728-6.
#' @examples
#' pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
#' sapply(pkgs, require, character.only = TRUE)
#' library(scMayoMap)
#' data(data)
#' obj <- scMayoMap(data = data, tissue = 'muscle')
#' plt <- scMayoMap.plot(scMayoMap.object = obj, directory = '~/Desktop/', width = 8, height = 6)

scMayoMap.plot <- function(scMayoMap.object, directory = NULL, width = 8, height = 6){
  annotation.mean <- scMayoMap.object$annotation.mean
  annotation.norm <- scMayoMap.object$annotation.norm
  top.n <- scMayoMap.object$res
  MaxCell.Cluster <- unique(as.vector(apply(annotation.mean, 1, function(x) names(unlist(sort(x, decreasing = T))[1:2]))))
  MaxCell.Cluster <- MaxCell.Cluster[!is.na(MaxCell.Cluster)]
  annotation.norm <- annotation.norm[,MaxCell.Cluster]
  annotation.norm$cluster <- rownames(annotation.norm)
  suppressMessages(mt.annotation <- melt(annotation.norm))
  colnames(mt.annotation) <- c('cluster','cell','score')
  if(!is.null(scMayoMap.object$tissue)){
    mt.annotation$cell <- gsub('.*:','',mt.annotation$cell)
    MaxCell.Cluster <- gsub('.*:','',MaxCell.Cluster)
  }

  ## prepare circle (top n data)
  tmp.df <- top.n[,c('celltype','ES.norm','cluster')] %>% dplyr::mutate(celltype = strsplit(celltype, "; ")) %>%
    tidyr::unnest(celltype) %>% dplyr::group_by(cluster) %>% dplyr::mutate(row = rev(dplyr::row_number())) %>% tidyr::spread(row, celltype)
  tmp.df1 <- tmp.df[,c('cluster','1')]
  colnames(tmp.df1)[2] <- 'cell'
  score1 <- mt.annotation %>% left_join(tmp.df1 %>% mutate(match = 't'), by = c('cluster','cell'))
  score1$score[is.na(score1$match)] <- NA
  score1 <- score1[,!(colnames(score1) %in% 'match')]

  if(ncol(tmp.df) <=3){
    cat('No biased celltype can be found.\n')
    plot.top = 1
  }else{
    tmp.df2 <- tmp.df[,c('cluster','2')]
    colnames(tmp.df2)[2] <- 'cell'
    score2 <- mt.annotation %>% left_join(tmp.df2 %>% mutate(match = 't'), by = c('cluster','cell'))
    score2$score[is.na(score2$match)] <- NA
    score2 <- score2[,!(colnames(score2) %in% 'match')]
    score2$cluster <- as.numeric(score2$cluster)
    plot.top = 2
  }

  mt.annotation$cluster <- as.numeric(gsub('Cluster ','', mt.annotation$cluster))
  score1$cluster <- as.numeric(gsub('Cluster ','', score1$cluster))
  ord <- unique((mt.annotation %>% group_by(cluster) %>% dplyr::slice_max(score, n=2))$cell)
  mt.annotation <- within(mt.annotation, cell <- factor(cell, levels= ord))
  score1 <- within(score1, cell <- factor(cell, levels= ord))
  suppressMessages(annotation.plot <- ggplot() +
                     geom_point(data = mt.annotation, aes(x=cluster,y=cell, size=score, color = score), shape =20) +
                     scale_color_continuous(type = "viridis") +
                     scale_x_continuous(breaks = seq(0, max(mt.annotation$cluster), by = 1)) +
                     geom_point(data = score1, aes(x=cluster,y=cell, size=score), shape = 21, color = 'red') +
                     theme_bw() +
                     labs(title= "", x="cluster", y="", caption = "Dot size also represents the score\nRed circle highlights top1 score\n Blue circle notates the 2nd top score") +
                     theme(axis.text = element_text(color = 'black', size = 14),
                           axis.text.x = element_text(color = 'black', size = 14),
                           axis.title = element_text(color = 'black', size = 14),
                           legend.title = element_text(color = 'black', size = 14),
                           legend.text = element_text(color = 'black', size = 14),
                           legend.position = 'right') +
                     guides(size = "none"))

  if(plot.top == 2){
    annotation.plot <- annotation.plot +geom_point(data = score2, aes(x=cluster,y=cell, size=score), shape = 21, color = 'blue')
  }

  suppressMessages(print(annotation.plot))
  if(is.null(directory)){directory <- getwd()}
  ggsave(paste0(directory, '/annotation.pdf'), width = width, height = height)

  return(annotation.plot)
}




