require(ggplot2)
require(dplyr)
require(tidyr)
require(tibble)
require(reshape2)

#' @title Calculate the cumulative variance of a vector
#' @param x a vector
#' @return {v}{a vector of cumulative variance}
#' @rdname cumvar
#' @import stats
#' @export
cumvar <- function(x){
  x <- x - x[sample.int(length(x), 1)]
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  return(v)
}


utils::globalVariables(c("celltype", "cluster", "cell","score"))


#' @title scMayoMap
#' @param data data frame. Rows: gene names; columns must include "avg_log2FC", "pct.1","pct.2", "gene".
#' @param database data frame. Already cleaned database saved in the package.
#' @param tissue tissue to be specified. If NULL, tissue will not be specified. User can choose one of:
#'                "adipose tissue"; "bladder"; "blood"; "bone"; "bone marrow";
#'                "brain";"breast"; "embryo"; "eye"; "gastrointestinal tract";
#'                "heart";"kidney";"liver";"lung";"mammary gland";"muscle";
#'                "other";"ovary";"pancreas";"placenta";"prostate";
#'                "skin";"spleen";"stomach";"testis";"thymus";"tooth";"uterus".
#' @param padj.cutoff number between 0-1. Used to filter the genes by FDR adjusted pvalue.
#' @param pct.cutoff number between 0-1. Genes with pct.1 greater than this number will be remained for later analysis.
#' @return A list with the elements
#' \item{res}{A data.frame listed cluster("cluster") with the top N potential celltypes("celltype") with estimated scores "ES.raw", standard errors "ES.se", normalized estimated scores "ES.norm", and markers matched for predicted cell type.}
#' \item{markers}{A data.frame of cell types scores of each cluster. Row: cluster, Column: annotated cell types.}
#' \item{tissue}{tissue specified in the analysis.}
#' \item{annotation.norm}{A data.frame of cell types scores of each cluster. Row: cluster, Column: annotated cell types.}
#' \item{annotation.mean}{A data.frame of cell types mean of each cluster. Row: cluster, Column: annotated cell types.}
#' @rdname scMayoMap
#' @importFrom dplyr summarise
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom tibble column_to_rownames
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @export scMayoMap
#'
scMayoMap <- function(data, database=NULL, padj.cutoff = 0.05, pct.cutoff = 0.25, tissue = NULL){
  if(is.null(database)){
    database = scMayoMapDatabase
  }

  ## Define tissue or not
  if(is.null(tissue)){
    db <- database[,!(colnames(database) %in% 'tissue')]
  }else{
    db <- (database[database$tissue==tissue,,drop =F])[,!(colnames(database) %in% 'tissue')]
  }
  data1 <- data
  # data1$cluster <- paste0('Cluster ', data1$cluster)
  colnames(data1)[grep('log',colnames(data1))] = "avg_log2FC"

  ## Filter gene by cutoffs
  data1 <- data1[data1$p_val_adj <= padj.cutoff & data1$pct.1 >= pct.cutoff,]
  data1$gene <- toupper(data1$gene)

  ## Match the data1base
  min.pct2 <- min(data1$pct.2[data1$pct.2>0])
  data1$Score <- ifelse(data1$pct.2 !=0, (2^data1$avg_log2FC * data1$pct.1)/data1$pct.2,(2^data1$avg_log2FC * data1$pct.1)/min.pct2)
  merged <- na.omit(merge(x=data1[,c('cluster','gene','Score')], y=db, all.y = TRUE))
  merged[,!(colnames(merged) %in% c('gene','Score','cluster'))] <- sapply(merged[,!(colnames(merged) %in% c('gene','Score','cluster'))], function(x) merged$Score * x )
  merged1 <- merged[,!(colnames(merged) %in% c('gene','Score')), drop =F]

  ## Calculate the mean value(dplyr is faster than baseR)
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

  ## save marker genes for mapped cell type for Seurat Dotplot
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




#' @title plot the annotation result of scMayoMap
#' @param scMayoMap.object return from function \code{scMayoMap}.
#' @param directory character. the directory to save the figures. Default is the current working directory.
#' @param width the width of the graphics region in inches. See R function \code{pdf}.
#' @param height the height of the graphics region in inches. See R function \code{pdf}.
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
#'

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
    cat('No biased celltype can be found.')
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
                     # labs(title= "", x="", y="", caption = "Dot size also represents the score\nRed circle highlights top1 score") +
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




