# Figure 3
# By Andrew Phan

rm(list=ls()) #clears variables
cat("\014") #clears console

library(CellChat)
library(Seurat)
library(patchwork)
library(future)
library(gplots)
library(viridis)
library(gplots)
library(reshape2)
library(ComplexHeatmap)
library(readr)
library(tidyverse)

#install.packages("Signac")

#load in Seurat object of combined lacZ and S3 scRNA-seq
load("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/Seurat/lacZ+S3_scRNAseq_reclustered.RData")
#LacZcellchat <- readRDS("/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/Changhan_Analysis/cellchat_LacZ.rds")
#S3cellchat <- readRDS("/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/Changhan_Analysis/cellchat_Schroom3.rds")
#load("/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/cellchat_finishedcomparativeanalysis_lacZvsS3.RData") #load in finished CellChat analysis
#load in finished CellChat analysis with all DBs and all timepoints combined. no split fibr.
load("~/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/ALLDB_L+S+combined.RDATA")

########### netAnalysis_signalingRole_heatmap_topbar ###############
#returns the top barplot matrix (each cell type incoming or outgoing signaling)
netAnalysis_signalingRole_heatmap_topbar <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                                     color.use = NULL, color.heatmap = "BuGn",
                                                     title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }

  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA


  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }

  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)

  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  draw(ht1)
  return(colSums(mat.ori))
}
##################

########### Split Fibr. Cell Type into Days to Run CellChat Analysis ###################
#input: a Seurat object. Should include metadata on Time and Cell Type (active identity)
#output: new CellChat object with identities as Cluster_Day (cell types, except Fibr. are split into D6/9/12/15) or specified group.by
seurat2CellChat_Fibr <- function(sObj, group.by = "Cluster_Day", splitFibr = T){
  data.input <- GetAssayData(object = sObj, assay = "RNA", slot = "data") # normalized data matrix
  meta = as.data.frame(sObj@meta.data) #extract current metadata
  meta$id = Idents(sObj) #get cell type labels

  if(splitFibr == T){ #if fibr. to be split by day, create new Cluster_Day metadata
    meta$Cluster_Day = 1:nrow(meta) #initialize Cluster_Day metadata column
    for(x in 1:nrow(meta)){ #create column of Cluster + Day in metadata
      if(meta$id[x]=="Fibr."){
        meta$Cluster_Day[x] = paste(meta$Time[x], meta$id[x], sep = "")
      }
      else{
        meta$Cluster_Day[x] = paste(meta$id[x])
      }
    }
    meta = meta[,-38] #get rid of incorrect cluster metadata
    sObj <- SetIdent(sObj, value = meta$Cluster_Day) #set the active identity to the new Cluster_Day column
  }
  else{ #else, just set the new metadata
    meta = meta[,-38] #get rid of incorrect cluster metadata
    sObj <- SetIdent(sObj, value = meta$id) #set the active identity to the new Cluster_Day column
  }
  show(UMAPPlot(sObj)) #plot UMAP to confirm new cell type groups

  cellchatreturn <- createCellChat(object = data.input, meta = meta, group.by = group.by)

  return(cellchatreturn)
}

#splitSeurat = SplitObject(combined2, split.by = "Condition")

#create LacZ CellChat object with Fibr. split into Days
#LacZcellchat <- seurat2CellChat_Fibr(splitSeurat$LacZ)
#create SHROOM3 CellChat object with Fibr. split into Days
#S3cellchat <- seurat2CellChat_Fibr(splitSeurat$Shroom3)

#All the calculations after `computeCommunProb` should be re-run!!
#  These include but not limited to `computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`.

########### Fixed Compute Similarity Related Functions ############
# source - https://github.com/sqjin/CellChat/issues/167

computeNetSimilarity_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  prob = methods::slot(object, slot.name)$prob
  if (is.null(k)) {
    if (dim(prob)[3] <= 25) {
      k <- ceiling(sqrt(dim(prob)[3]))
    } else {
      k <- ceiling(sqrt(dim(prob)[3])) + 1
    }
  }
  if (!is.null(thresh)) {
    prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
  }
  if (type == "functional") {
    # compute the functional similarity
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    S2 <- D_signalings; S3 <- D_signalings;
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi + Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
    # S_signalings <- S1 *S2
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    D_signalings <- D_signalings + t(D_signalings)
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- dimnames(prob)[[3]]
  colnames(Similarity) <- dimnames(prob)[[3]]
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}


computeNetSimilarityPairwise_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, thresh = NULL) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("Compute signaling network similarity for datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  net <- list()
  signalingAll <- c()
  object.net.nameAll <- c()
  # 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    object.net.name <- names(methods::slot(object, slot.name))[comparison[i]]
    object.net.nameAll <- c(object.net.nameAll, object.net.name)
    net[[i]] = object.net$prob
    signalingAll <- c(signalingAll, paste0(dimnames(net[[i]])[[3]], "--", object.net.name))
    # signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
  }
  names(net) <- object.net.nameAll
  net.dim <- sapply(net, dim)[3,]
  nnet <- sum(net.dim)
  position <- cumsum(net.dim); position <- c(0,position)
  if (is.null(k)) {
    if (nnet <= 25) {
      k <- ceiling(sqrt(nnet))
    } else {
      k <- ceiling(sqrt(nnet)) + 1
    }
  }
  if (!is.null(thresh)) {
    for (i in 1:length(net)) {
      neti <- net[[i]]
      neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
      net[[i]] <- neti
    }
  }
  if (type == "functional") {
    # compute the functional similarity
    S3 <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi + Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0;  diag(S3) <- 1
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- signalingAll
  colnames(Similarity) <- rownames(Similarity)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  # methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}

### netEmbedding
netEmbedding_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, pathway.remove = NULL, k = NULL) {
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Manifold learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Manifold learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
  if (is.null(pathway.remove)) {
    pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
  }
  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
    Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
  }
  if (is.null(k)) {
    k <- ceiling(sqrt(dim(Similarity)[1])) + 1
  }
  options(warn = -1)
  # dimension reduction
  Y <- runUMAP(Similarity)#, min.dist = 0.3, n.neighbors = k) #threw an error
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$dr)) {
    methods::slot(object, slot.name)$similarity[[type]]$dr <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
  return(object)
}



netClustering_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, methods = "kmeans", do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, k.eigen = NULL) {
  type <- match.arg(type)
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Classification learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Classification learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Y[is.na(Y)] <- 0
  data.use <- Y
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use,k,nstart=10)$cluster
    } else {
      N <- nrow(data.use)
      kRange <- seq(2,min(N-1, 10),by = 1)
      if (do.parallel) {
        future::plan("multiprocess", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
      )
      results = my.sapply(
        X = 1:length(kRange),
        FUN = function(x) {
          idents <- kmeans(data.use,kRange[x],nstart=10)$cluster
          clusIndex <- idents
          #adjMat0 <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
          adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
          return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
        },
        simplify = FALSE
      )
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce('+', adjMat)/length(kRange)
      res <- computeEigengap(as.matrix(CM))
      numCluster <- res$upper_bound
      clusters = kmeans(data.use,numCluster,nstart=10)$cluster
      if (do.plot) {
        gg <- res$gg.obj
        ggsave(filename= paste0("estimationNumCluster_",fig.id,"_",type,"_dataset_",comparison.name,".pdf"), plot=gg, width = 3.5, height = 3, units = 'in', dpi = 300)
      }
    }
  } else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D)-A                       # unnormalized version
    L <- diag(D^-0.5)%*%L%*% diag(D^-0.5) # normalized version
    evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
    # pick the first k first k eigenvectors (corresponding k smallest) as data points in spectral space
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[,(ncol(evL$vectors)-k.eigen1):ncol(evL$vectors)]
    clusters = kmeans(Z,k,nstart=20)$cluster
  }
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
    methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
  return(object)
}
#############

######## Changhan's Network Similarity Barplot Function ########
SimilarityCompute <- function (object, slot.name = "netP", type = c("functional", "structural"), comparison1 = NULL,
                               comparison2 = c(1,2), x.rotation = 90, title = NULL, color.use = NULL, bar.w = NULL, font.size = 8){
  type <- match.arg(type)
  if (is.null(comparison1)) {
    comparison1 <- 1:length(unique(object@meta$datasets))
  }
  comparison.name <- paste(comparison1, collapse = "-")
  cat("Compute the distance of signaling networks between datasets",
      as.character(comparison1[comparison2]), "\n")
  comparison2.name <- names(methods::slot(object, slot.name))[comparison1[comparison2]]
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  data1 <- Y[group %in% comparison2.name[1], ]
  data2 <- Y[group %in% comparison2.name[2], ]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))
  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))
  dist <- NULL
  for (i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i, ],
                                               data2[i, ])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))
  return(df)

}
################

object.list <- list(L = LacZcellchat, S = S3cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

############ INDIVIDUAL CELLCHAT ANALYSIS ####################
########### Convert Seurat object to CellChat object ############

library(Seurat) #need Seurat package loaded for SplitObject and GetAssayData

splitSeurat = SplitObject(combined2, split.by = "Condition")

#create LacZ CellChat object
data.input <- GetAssayData(object = splitSeurat$LacZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(splitSeurat$LacZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#LacZcellchat <- createCellChat(object = splitSeurat$LacZ, assay = "RNA", meta = meta, group.by = "group")
LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#create SHROOM3 CellChat object
data.input <- GetAssayData(object = splitSeurat$Shroom3, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(splitSeurat$Shroom3)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

S3cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#if running on other cellchat objects
LacZcellchat = D9LacZcellchat
S3cellchat = D6S3cellchat
LacZcellchat@idents <- droplevels(LacZcellchat@idents, exclude = setdiff(levels(meta$labels),unique(meta$labels))) #if idents labels are not all used, need to drop unused
S3cellchat@idents <- droplevels(S3cellchat@idents, exclude = setdiff(levels(meta$labels),unique(meta$labels)))

######### Set the ligand-receptor interaction database ############

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

#CellChatDB.use <- subsetDB(CellChatDB, search = c("ECM-Receptor", "Cell-Cell Contact")) #focus on cell-ECM and cell-cell interactions (not para/autocrine)

CellChatDB.use <- subsetDB(CellChatDB, search = c("ECM-Receptor", "Cell-Cell Contact", "Secreted Signaling")) #use all signaling
LacZcellchat@DB <- CellChatDB.use
S3cellchat@DB <- CellChatDB.use

############ Preprocessing the expression data for cell-cell communication analysis ############

LacZcellchat <- subsetData(LacZcellchat) # This step is necessary even if using the whole database
S3cellchat <- subsetData(S3cellchat) # This step is necessary even if using the whole database

future::plan("multiprocess", workers = 4) # do parallel

LacZcellchat <- identifyOverExpressedGenes(LacZcellchat)
LacZcellchat <- identifyOverExpressedInteractions(LacZcellchat)

S3cellchat <- identifyOverExpressedGenes(S3cellchat)
S3cellchat <- identifyOverExpressedInteractions(S3cellchat)

#cellchat <- projectData(cellchat, PPI.human) # project gene expression data onto PPI network (optional)

############ Inference of cell-cell communication network ##########

LacZcellchat <- computeCommunProb(LacZcellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
LacZcellchat <- filterCommunication(LacZcellchat, min.cells = 10)

S3cellchat <- computeCommunProb(S3cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
S3cellchat <- filterCommunication(S3cellchat, min.cells = 10)

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

#Infer the cell-cell communication at a signaling pathway level
LacZcellchat <- computeCommunProbPathway(LacZcellchat)
S3cellchat <- computeCommunProbPathway(S3cellchat)

#Calculate the aggregated cell-cell communication network
LacZcellchat <- aggregateNet(LacZcellchat)
S3cellchat <- aggregateNet(S3cellchat)

#visualize the aggregated cell-cell communication network by showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
# groupSize <- as.numeric(table(LacZcellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(LacZcellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(LacZcellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# examine the signaling sent from each cell group. also control the parameter edge.weight.max so that we can compare edge weights between different networks
# mat <- LacZcellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

################### Visualization of cell-cell communication network ####################


#see signaling pathways with significant interaction
LacZcellchat@netP$pathways
S3cellchat@netP$pathways

# pathways.show <- c("FN1")
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
# vertex.receiver = seq(1,4) # a numeric vector.
# netVisual_aggregate(LacZcellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
#
# # Circle plot
# par(mfrow=c(1,1))
# netVisual_aggregate(LacZcellchat, signaling = pathways.show, layout = "circle")
#
# # Chord diagram
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#
# # Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# #> Do heatmap based on a single object
#
# # Chord diagram
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
# #> Plot the aggregated cell-cell communication network at the signaling pathway level

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
# netAnalysis_contribution(LacZcellchat, signaling = pathways.show)

# #SEE ALL GRAPHS FOR PATHWAYS
# # Access all the signaling pathways showing significant communications
# pathways.show.all <- cellchat@netP$pathways
# # check the order of cell identity to set suitable vertex.receiver
# levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

################## Systems analysis of cell-cell communication network ###################

# Compute the network centrality scores
LacZcellchat <- netAnalysis_computeCentrality(LacZcellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#netAnalysis_signalingRole_network(LacZcellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Compute the network centrality scores
S3cellchat <- netAnalysis_computeCentrality(S3cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#netAnalysis_signalingRole_network(S3cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space

# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# gg1 <- netAnalysis_signalingRole_scatter(LacZcellchat)
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2 <- netAnalysis_signalingRole_scatter(LacZcellchat, signaling = c("FN1"))
# #> Signaling role analysis on the cell-cell communication network from user's input
# gg1 + gg2

# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# gg1 <- netAnalysis_signalingRole_scatter(S3cellchat)
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2 <- netAnalysis_signalingRole_scatter(S3cellchat, signaling = c("FN1"))
# #> Signaling role analysis on the cell-cell communication network from user's input
# gg1 + gg2

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NCCs", signaling.exclude = c("COLLAGEN", "LAMININ", "CDH", "MK"))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "iPSCs", signaling.exclude = c("COLLAGEN", "LAMININ", "CDH", "MK"))
gg1 + gg2

# #Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "incoming")
ht1 + ht2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(S3cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(S3cellchat, pattern = "incoming")
ht1 + ht2

# #Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# ht1 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "outgoing")
# ht2 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "incoming")
# ht1 + ht2

#UNCOMMENT THE BELOW SECTION LATER
#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
# library(NMF)
# library(ggalluvial)
#
# #NOTE: NEED TO ADJUST nPatterns PER CONDITION
# selectK(LacZcellchat, pattern = "outgoing")
# nPatterns = 5 #change to before where both graphs drop off together. ECM+cell-cell
# nPatterns = 3 #all DB
# LacZcellchat <- identifyCommunicationPatterns(LacZcellchat, pattern = "outgoing", k = nPatterns)
# netAnalysis_river(LacZcellchat, pattern = "outgoing")
#
# selectK(S3cellchat, pattern = "outgoing")
# nPatterns = 2
# nPatterns = 3 #all DB
# S3cellchat <- identifyCommunicationPatterns(S3cellchat, pattern = "outgoing", k = nPatterns)
# netAnalysis_river(S3cellchat, pattern = "outgoing")
#
# selectK(LacZcellchat, pattern = "incoming")
# nPatterns = 3
# nPatterns = 2 #all DB
# LacZcellchat <- identifyCommunicationPatterns(LacZcellchat, pattern = "incoming", k = nPatterns)
# netAnalysis_river(LacZcellchat, pattern = "incoming")
#
# selectK(S3cellchat, pattern = "incoming")
# nPatterns = 2 #was 3 before
# nPatterns = 3 #all DB
# S3cellchat <- identifyCommunicationPatterns(S3cellchat, pattern = "incoming", k = nPatterns)
# netAnalysis_river(S3cellchat, pattern = "incoming")

#Manifold and classification learning analysis of signaling networks
LacZcellchat <- computeNetSimilarity_GitHub(LacZcellchat, type = "functional")
LacZcellchat <- netEmbedding_GitHub(LacZcellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
LacZcellchat <- netClustering_GitHub(LacZcellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(LacZcellchat, type = "functional", label.size = 3.5)

LacZcellchat <- computeNetSimilarity_GitHub(LacZcellchat, type = "structural")
LacZcellchat <- netEmbedding_GitHub(LacZcellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
LacZcellchat <- netClustering_GitHub(LacZcellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(LacZcellchat, type = "structural", label.size = 3.5)

# Y = S3cellchat
#
# Y <- methods::slot(S3cellchat, data)$similarity[[type]]$dr[[comparison.name]]
# Y <- Y[-is.nan(Y[,1]), ]
# methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y

S3cellchat <- computeNetSimilarity_GitHub(S3cellchat, type = "functional")
S3cellchat <- netEmbedding_GitHub(S3cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
S3cellchat <- netClustering_GitHub(S3cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(S3cellchat, type = "functional", label.size = 3.5)

S3cellchat <- computeNetSimilarity(S3cellchat, type = "structural")
S3cellchat <- netEmbedding(S3cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
S3cellchat <- netClustering(S3cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(S3cellchat, type = "structural", label.size = 3.5)

#saveRDS(LacZcellchat, file = "/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/CellChat/LacZcellchat_individual.rds")
#saveRDS(S3cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/CellChat/S3cellchat_individual.rds")

#saveRDS(LacZcellchat, file = "/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_lacZ.RDS")
#saveRDS(S3cellchat, file = "/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_S3.RDS")

##################### Merge the LacZ and S3 CellChat Objects After Individual Analysis #######################
#create a merged CellChat object of LacZ and S3 conditions
object.list <- list(lacZ = LacZcellchat, S3 = S3cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#split fibr. all days combined
object.list <- list(lacZ = ALLDB_lacZ_splitFibr, S3 = ALLDB_S3_splitFibr)
combined <- mergeCellChat(object.list, add.names = names(object.list))

#all databases
object.list <- list(lacZ = D6ALLDB_lacZ, S3 = D6ALLDB_S3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#all databases liftover
object.list <- list(lacZ = D6ALLDB_lacZ_LIFTUP, S3 = D6S3cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat = cellchatD6
cellchat = cellchatD9
cellchat = cellchatD12
cellchat = cellchatD15

#saveRDS(combined, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/ALLDB_L+Scellchat_combined_splitFibr.rds")

############### COMPARATIVE CELLCHAT ANALYSIS ##################################################################
############### Predict general principles of cell-cell communication ###############
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)#, sources.use = c(2:5,7), targets.use = c(2:5,7), remove.isolate = T)


#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight") #cannot specify specific cell types: sources.use = c(2:5,7), targets.use = c(2:5,7), remove.isolate = T
#red/blue colored edges represent increased/decreased signaling in the 2nd dataset compared to the 1st

gg1 <- netVisual_heatmap(cellchat)#, row.show = c(1,8,11,4,5,2,3,6,7,13,10,12,9), col.show = c(1,8,11,4,5,2,3,6,7,13,10,12,9))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")#, row.show = c(1,8,11,4,5,2,3,6,7,13,10,12,9), col.show = c(1,8,11,4,5,2,3,6,7,13,10,12,9))
#> Do heatmap based on a merged object
gg1 + gg2


#To better control the node size and edge weights of the inferred networks across different datasets,
#we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Differential number of interactions or interaction strength among different cell types
#To simplify the complicated network and gain insights into the cell-cell communication at the cell type level
#we can aggregate the cell-cell communication based on the defined groups of cells (ex: prolif. troph. + troph. = trophs.)
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
# group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
# object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

############## Identify the conserved and context-specific signaling pathways ###############
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
cellchat <- computeNetSimilarityPairwise_GitHub(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding_GitHub(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering_GitHub(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#Error in data.frame(x = Y[, 1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),  :
#                      arguments imply differing number of rows: 60, 0

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise_GitHub(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding_GitHub(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering_GitHub(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")

#Extract pathway distances per signaling pathway and plot heatmap
rankPlot <- rankSimilarity(cellchat, type = "functional")
rank = data.matrix(rankPlot$data)
rank[,1] = 1
rank = log(x=rank, base = 10)
#plot log(distance) of pathways
heatmap.2(t(rank), Colv = NA, Rowv = NA, col = magma(256))

#Identify and visualize the conserved and context-specific signaling pathways
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)

# combining all the identified signaling pathways from different datasets
i=1
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
pathway.union = c("FN1", "ncWNT", "COLLAGEN")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

################### Identify the upregulated and down-regulated signaling ligand-receptor pairs ######################
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#FN1 pairs downregulated communication in S3 interactions?

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LacZ", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LacZ", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LacZ"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LacZ
net.up <- subsetCommunication(cellchat, net = net, datasets = "LacZ",ligand.logFC = 0.1, receptor.logFC = NULL)
#NOTE: Optimize ligand.logFC and receptor.logFC

#Error in subsetCommunication_internal(net, LR, cells.level, slot.name = slot.name,  :
#   No significant signaling interactions are inferred based on the input!

# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in S3, i.e.,downregulated in LacZ
net.down <- subsetCommunication(cellchat, net = net, datasets = "S3",ligand.logFC = -0.1, receptor.logFC = -0.1)

#see individual signaling genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:11), targets.use = c(1:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:11), targets.use = c(1:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#Error in netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11),  :
#   No signaling links are inferred!

############### Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram ####################
pathways.show <- c("FN1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("ncWNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], sources.use = c(7), targets.use = c(2:5,7), remove.isolate = T, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], signaling.name = paste(pathways.show, names(object.list)[i]))
}

#draw individual heatmaps for LacZ and S3 condition ncWNT signaling
pathways.show <- c("ncWNT")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), color.heatmap = c("white", "darkred"))
pathways.show=""
netVisual_heatmap(combined, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), color.heatmap = c("white", "darkred"), sources.use = c(3,4,7), targets.use = c(3,4,7), remove.isolate = T)

#draw heatmap of ncWNT comparison
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(2:5,7), targets.use = c(1:10))
#heatmap of NCCs, SMCs, and iPSCs. measure = weight
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(2,4,6,9), targets.use = c(2,4,6,9), remove.isolate = T, cluster.rows = F)
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(1:9), targets.use = c(1:9), remove.isolate = T, cluster.rows = F)
ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))

netVisual_heatmap(cellchat, comparison = c(1,2), signaling = "FN1", measure = "weight", title.name = paste(pathways.show, "signaling "), color.heatmap = c("white", "darkred"))

mat = ht@matrix

#custom heatmap to adjust limits
library(circlize)
library(RColorBrewer)
col_fun = colorRamp2(c(-.9, 0, .9), c("darkblue", "white", "red")) #create a linear color palette from the min and max of the colors
Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F, heatmap_legend_param = list(at = pretty(c(-1, 1)), color_bar = "continuous"), col = col_fun)#viridis(option = "", n = 10))

#Split Fibr. Analyses:
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(1:13), targets.use = c(1:13))
ht <- netVisual_heatmap(cellchat, color.heatmap = c("white", "darkred"), signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(2:5,7), targets.use = c(2:5,7), remove.isolate = T, row.show = c(3,4,1,2,5), col.show = c(3,4,1,2,5))
ht <- netVisual_heatmap(cellchat, color.heatmap = c("white", "darkred"), signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(3,4,7), targets.use = c(3,4,7), remove.isolate = T, row.show = c(2,1,3), col.show = c(2,1,3))

ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))

col_fun = colorRamp2(c(-.5, 0, 1.5), c("darkblue", "white", "red")) #create a linear color palette from the min and max of the colors
Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F, heatmap_legend_param = list(at = pretty(c(-1, 1)), color_bar = "continuous"), col = col_fun)#viridis(option = "", n = 10))

col_fun = colorRamp2(c(0, 1.3), c("white", "darkred")) #create a linear color palette from the min and max of the colors
Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F, heatmap_legend_param = list(at = pretty(c(0, 1.5)), color_bar = "continuous"), col = col_fun)#viridis(option = "", n = 10))


#draw heatmap of FN1 comparison
pathways.show <- c("WNT")
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(2:5,7), targets.use = c(1:10))
#heatmap of NCCs, SMCs, and iPSCs. measure = weight
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(2,4,6,9), targets.use = c(2,4,6,9), remove.isolate = T, cluster.rows = F)
ht <- netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight", title.name = paste(pathways.show, "signaling "), sources.use = c(1:10), targets.use = c(1:10), remove.isolate = F, cluster.rows = F)
ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))


# Chord diagram
# pathways.show <- c("ncWNT")
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
# }

# # Aggregated Chord diagram
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(object.list[[1]]@idents)
# pathways.show <- c("CXCL")
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
# }
# # Plot the aggregated cell-cell communication network at the signaling pathway level

#Chord diagram of single signaling pathway for given sources/targets
# par(mfrow = c(1, 2), xpd=TRUE)
# # compare all the interactions sending from subpops to iPSCs or NCCs cells
# for (i in 1:length(object.list)) {
#   netVisual_chord_cell(object.list[[1]], signaling = "NOTCH", sources.use = c(2, 3, 5, 6), targets.use = c(2, 3, 5, 6), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
# }

#`0` = "CIMHC+", `1` = "NCCs", `2` = "Prolif.Troph.", `3` = "Fibr.", `4` = "Epi.", `5` = "iPSCs",  `6` = "Troph.", `7` = "Neur.Epi.", `8` = "SMCs", `9` = "Neur."

#FN1 signaling reduced in NCCs, iPSCs, SMCs, Epithelial, Prolif. Troph.
pathways.show <- c("FN1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(2,3,5,6,9), targets.use = c(2,3,5,6,9), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#EPHA signaling increased in NCCs, iPSCs, SMCs, Epithelial, Prolif. Troph.
pathways.show <- c("EPHA")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(2,6,9), targets.use = c(2,6,9), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#NOTCH signaling reduced in NCCs, iPSCs, SMCs, Epithelial, Prolif. Troph.
pathways.show <- c("NOTCH")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(2,3,5,6,9), targets.use = c(2,3,5,6,9), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#CIRCLE PLOTS
#ncWNT outgoing
pathways.show <- c("ncWNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, label.edge = T, remove.isolate = T, sources.use = c(2,6,9), targets.use = c(2,6,9), layout = "circle", edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#ncWNT incoming
pathways.show <- c("ncWNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(1:9), targets.use = c(2), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1, 2), xpd=TRUE)
#1 = CIMHC+ 2 = NCC 3 = prolif. troph. 4 = fibr. 5 = epi. 6 = iPSCs 7 = troph. 8 = neur. epi. 9 = SMCs 10 = Neur.
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], signaling = "LAMININ", sources.use = c(1), targets.use = c(1), slot.name = "netP", title.name = paste0("Sig", names(object.list)[i]), legend.pos.x = 10)
}

#for split fibr.
pathways.show <- c("ncWNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(2:5,7), targets.use = c(2:5,7), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

group.cellType <- c(2,1,1,2,1,1,2,1,2,1) # grouping cell clusters based on S3 in
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("FN1")
par(mfrow = c(1,1), xpd=TRUE, mar = c(1, 1, 1, 1))
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], remove.isolate = T, sources.use = c(1,2,3,4,5,6,8,9,10), targets.use = c(1,2,3,4,5,6,8,9,10), signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

#Vertex FIGURES:
#visualize vertex graph with iPSCs, NCCs, Prolif. Troph., and Epi. cells on left
vertex.receiver = c(2,3,5,6,9)
pathways.show = "NOTCH"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

vertex.receiver = c(2,3,5,6)
pathways.show = "ncWNT"
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show))

# Aggregated Chord diagram
#LacZ in
# 1: 9
# 2: 2, 3, 5, 6
# 3: 1, 4, 8, 10
group.cellType <- c(3,2,2,3,2,2,0,3,1,3) # grouping cell clusters based on LacZ in
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("LAMININ")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

#S3 in
# 1: 2, 3, 5, 6, 8, 10
# 2: 1, 4, 7, 9
group.cellType <- c(2,1,1,2,1,1,2,1,2,1) # grouping cell clusters based on S3 in
names(group.cellType) <- levels(object.list[[1]]@seurat.clusters)
pathways.show <- c("NOTCH")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

#LacZ out
#1: 4, 9
#2: 1, 3, 5, 7, 10
#3: 6
#4: 8
group.cellType <- c(2,0,2,1,2,3,2,4,1,2) # grouping cell clusters based on LacZ out
names(group.cellType) <- levels(object.list[[1]]@seurat.clusters)
pathways.show <- c("NOTCH")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

#S3 out
#1: 1, 4, 8, 9
#2: 2, 3, 5, 6, 7, 10
group.cellType <- c(1,2,2,1,2,2,2,1,1,2) # grouping cell clusters based on S3 out
names(group.cellType) <- levels(object.list[[1]]@seurat.clusters)
pathways.show <- c("NOTCH")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

pathways.show <- c("ncWNT")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

netVisual_diffInteraction(object = cellchat, measure = "weight", sources.use = c(2:5,7), targets.use = c(2:5,7), remove.isolate = T, )
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)

#ncWNT interactions:
#heatmap of senders and receivers
pathways.show <- c("ncWNT")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ht[[i]] <- netVisual_heatmap(cellchat, signaling = pathways.show, title.name = paste(pathways.show, "signaling "), measure = "weight")
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

################# Compare the signaling gene expression distribution between different datasets ################
# cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("L", "S")) # set factor level
plotGeneExpression(cellchat, signaling = "ncWNT", split.by = "Cluster_Day", colors.ggplot = T, split.plot = T)
plotGeneExpression(cellchat, signaling = "NOTCH", colors.ggplot = T, split.plot = F)

StackedVlnPlot(object = cellchat, features = "EPHA")

#saveRDS(cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/CellChat/cellchat_finishedcomparisonAnalysis_lacZvsS3.rds")

############## Analyze Day by Day Cell Signaling Interactions ###############

# cell.use = unlist(cellchat@idents)[cellchat@meta$Time == "D12"]
# test = subsetCellChat(cellchat, cells.use = cell.use)
#
# cellchat2 <- setIdent(cellchat2, ident.use = "Time")
# test = subsetCellChat(cellchat2, idents.use = c("D12", "D15"))
#
# cellchat2 = cellchat
# test = subsetCellChat(cellchat2, idents.use = c(1))

# D12cellchat = subset(x = cellchat, subset = meta$Time == "D12")

i=1
i=2
vertex.receiver = c(2,3,5,6,9)
pathways.show = "FN1"
for (i in 1:2) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#Non-split Fibr.
object.list <- list(lacZ = LacZcellchat, S3 = S3cellchat)
pathways.show = "ncWNT"
weight.max <- getMaxWeight(object.list[1:2], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfcol = c(1,2), xpd=TRUE)
#groups
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           targets.use = c(2,4,6,9), sources.use = c(2,4,6,9), remove.isolate = T, top = 1, label.edge=T))
}

#all
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           top = .25))
}

#individual
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           sources.use = c(4), targets.use = c(1:9), remove.isolate = T, top = 1, label.edge=T))
}

dev.off()
dev.new()

object.list <- list(lacZ = ALLDB_lacZ_splitFibr, S3 = ALLDB_S3_splitFibr)
#Split Fibr.
pathways.show = "ncWNT"
weight.max <- getMaxWeight(object.list[1:2], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfcol = c(1,2), xpd=TRUE)
#groups
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           targets.use = c(2:5,7,8,12), sources.use = c(2:5,7,8,12), remove.isolate = T)) #can add top = 0.1 and label.edge = T to look at top 10% and add labels to strength
}

#all
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           top = .25))
}

#individual
#2=D12 fibr, 3=D15, 4=D6, 5=D9, 6=epi, 7=iPSCs, 8=NCCs, 9=neur. 10=neur.epi., 11=prolif.troph., 12=SMCs
for (i in 1:2) {
  show(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),alpha.edge = 0.8,
                           sources.use = c(1:12), targets.use = c(12), remove.isolate = T))
}

dev.off()
dev.new()

#info flow stacked bar chart for NCCs + SMCs. FN1 and ncWNT both are enriched in LacZ
rankNet(cellchat, sources.use = c(1,9), targets.use = c(1,9), stacked = T)

############ Create CellChat Objects per Day ###############
# D6cellchat <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/CellChat_Day_Objs_cell-cell+ECMonly/cellchat_D6.rds")
# D9cellchat <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/CellChat_Day_Objs_cell-cell+ECMonly/cellchat_D9.rds")
# D12cellchat <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/CellChat_Day_Objs_cell-cell+ECMonly/cellchat_D12.rds")
# D15cellchat <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/CellChat_Day_Objs_cell-cell+ECMonly/cellchat_D15.rds")

#Create CellChat objects per day:
#split the seurat by day
combined2 = AddMetaData(combined2, combined2@active.ident, col.name = "celltype") #add active idents to metadata
splitSeurat = SplitObject(combined2, split.by = "Time")
condsplitSeurat = SplitObject(splitSeurat$D15, split.by = "Condition")

#test =
subset(x = splitSeurat$D6, cell.use = WhichCells(splitSeurat$D6, ident = "CIMHC+"))

UMAPPlot(splitSeurat$D9)
UMAPPlot(condsplitSeurat$LacZ)
UMAPPlot(condsplitSeurat$Shroom3)

subset(x = splitSeurat$D6, subset = (celltype == "CIMHC+" | celltype == "Fibr." | celltype == "Troph." | celltype == "Neur.Epi.")) #when subsetting with Seurat, make sure to separate subsets, or it will give you much fewer than expected. See https://www.biostars.org/p/420831/
subset(x = splitSeurat$D9, subset = (celltype == "CIMHC+" | celltype == "NCCs" | celltype == "Prolif.Troph." | celltype == "Fibr." | celltype == "Epi." | celltype == "Troph." | celltype == "Neur.Epi."| celltype == "SMCs")) #when subsetting with Seurat, make sure to separate subsets, or it will give you much fewer than expected. See https://www.biostars.org/p/420831/

ncol(x = test)
ncol(splitSeurat$D6)

#create the combined cellchat objs. need to split by condition and run the individual analyses separately before combining and running combined analysis

#need to subset D6 and D9 to only shared cell types. D12 and D15 have the same cell types
seuratD6 = subset(x = splitSeurat$D6, subset = (celltype == "CIMHC+" | celltype == "Fibr." | celltype == "Troph." | celltype == "Neur.Epi.")) #when subsetting with Seurat, make sure to separate subsets, or it will give you much fewer than expected. See https://www.biostars.org/p/420831/
seuratD9 = subset(x = splitSeurat$D9, subset = (celltype == "CIMHC+" | celltype == "NCCs" | celltype == "Prolif.Troph." | celltype == "Fibr." | celltype == "Epi." | celltype == "Troph." | celltype == "Neur.Epi."| celltype == "SMCs")) #when subsetting with Seurat, make sure to separate subsets, or it will give you much fewer than expected. See https://www.biostars.org/p/420831/

#without lift up
cellchatD6 = seurat2CellChat_Fibr(seuratD6, group.by = "celltype", splitFibr = F)
cellchatD9 = seurat2CellChat_Fibr(seuratD9, group.by = "celltype", splitFibr = F)
cellchatD12 = seurat2CellChat_Fibr(splitSeurat$D12, group.by = "celltype", splitFibr = F)
cellchatD15 = seurat2CellChat_Fibr(splitSeurat$D15, group.by = "celltype", splitFibr = F)

UMAPPlot(condsplitSeurat$Shroom3)

#for individual analyses: split the day cellchat objects into individual LacZ and SHROOM3 conditions
data.input = cellchatD6@data
meta = cellchatD6@meta
id = cellchatD6@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "LacZ"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D6LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD6@data
meta = cellchatD6@meta
id = cellchatD6@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "Shroom3"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D6S3cellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD9@data
meta = cellchatD9@meta
id = cellchatD9@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "LacZ"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D9LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD9@data
meta = cellchatD9@meta
id = cellchatD9@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "Shroom3"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D9S3cellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD12@data
meta = cellchatD12@meta
id = cellchatD12@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "LacZ"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D12LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD12@data
meta = cellchatD12@meta
id = cellchatD12@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "Shroom3"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D12S3cellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD15@data
meta = cellchatD15@meta
id = cellchatD15@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "LacZ"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D15LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

data.input = cellchatD15@data
meta = cellchatD15@meta
id = cellchatD15@idents
meta = cbind(meta, id)
cell.use = rownames(meta)[meta$Condition == "Shroom3"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
D15S3cellchat <- createCellChat(object = data.input, meta = meta, group.by = "id")

#NEXT RUN INDIVIDUAL ANALYSES, COMBINE AND THEN RUN COMPARATIVE ANALYSES. THEN CONTINUE TO NEXT SECTION TO FIND NORMALIZED FUNCTIONAL SIMILARITY VALUES
##################

#Find Normalized Functional Similarity Values

#load in individual cellchat objects that have been processed individually
D6ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D6ALLDB_lacZ.RDS")
D6ALLDB_S3 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D6ALLDB_S3.RDS")

D9ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_lacZ.RDS")
D9ALLDB_S3 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_S3.RDS")

D12ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D12ALLDB_lacZ.RDS")
D12ALLDB_S3 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D12ALLDB_S3.RDS")

D15ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D15ALLDB_lacZ.RDS")
D15ALLDB_S3 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D15ALLDB_S3.RDS")

#load in cellchat objects that contain LacZ and S3 conditions combined per timepoint
cellchatD6 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D6ALLDB_COMBINED.RDS")
cellchatD9 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_COMBINED.RDS")
cellchatD12 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D12ALLDB_COMBINED.RDS")
cellchatD15 <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D15ALLDB_COMBINED.RDS")

d6similarity = SimilarityCompute(cellchatD6, type = "functional")
d9similarity = SimilarityCompute(cellchatD9, type = "functional")
d12similarity = SimilarityCompute(cellchatD12, type = "functional")
d15similarity = SimilarityCompute(cellchatD15, type = "functional")

pathways.show = c("FN1", "ncWNT", "EPHA", "NOTCH")

similarity_matrix = data.frame(matrix(nrow = 4, ncol = length(pathways.show)+1))
colnames(similarity_matrix) = c("Time",pathways.show)
similarity_matrix$Time = c(6,9,12,15)

similarity_matrix[1,2:ncol(similarity_matrix)] = d6similarity[pathways.show,2]
similarity_matrix[2,2:ncol(similarity_matrix)] = d9similarity[pathways.show,2]
similarity_matrix[3,2:ncol(similarity_matrix)] = d12similarity[pathways.show,2]
similarity_matrix[4,2:ncol(similarity_matrix)] = d15similarity[pathways.show,2]

norm_sim_mat = similarity_matrix/max(similarity_matrix[,2:ncol(similarity_matrix)]) #normalize to max value of signaling pathways
norm_sim_mat$Time = similarity_matrix$Time #revert time to proper days again

# # Specify id.vars: the variables to keep but not split apart on
# df = melt(norm_sim_mat, id.vars=c("Time"))

ggplot(data=melt(norm_sim_mat, id.vars=c("Time")), aes(x=Time, y=value, fill=variable)) +
  geom_bar(stat="identity",  position=position_dodge()) + theme_classic()

#revert name of processed cellchat object
cellchatD6 = cellchat
cellchatD9 = cellchat
cellchatD12 = cellchat
cellchatD15 = cellchat

saveRDS(cellchatD9, file = "/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D12ALLDB_COMBINED_LIFTOVER.RDS")

#Lift Up for Comparative
#create combined cellchat objects
cellchatD6 = seurat2CellChat_Fibr(splitSeurat$D6, group.by = "celltype", splitFibr = F)
cellchatD9 = seurat2CellChat_Fibr(splitSeurat$D9, group.by = "celltype", splitFibr = F)
cellchatD12 = seurat2CellChat_Fibr(splitSeurat$D12, group.by = "celltype", splitFibr = F)
cellchatD15 = seurat2CellChat_Fibr(splitSeurat$D15, group.by = "celltype", splitFibr = F)
#create individual condition cellchat objects above

#lift up S3 conditions using LacZ idents
group.new = levels(D6ALLDB_lacZ_LIFTUP@idents)
D6S3cellchat <- liftCellChat(D6ALLDB_S3, group.new)

group.new = levels(D9ALLDB_lacZ_LIFTOVER@idents)
D9S3cellchat <- liftCellChat(D9ALLDB_S3, group.new)

d6similarity = SimilarityCompute(`D6ALLDB_COMBINED_ LIFT`, type = "functional")
d9similarity = SimilarityCompute(D9ALLDB_COMBINED_LIFTOVER, type = "functional")
d12similarity = SimilarityCompute(D12ALLDB_COMBINED, type = "functional")
d15similarity = SimilarityCompute(D15ALLDB_COMBINED, type = "functional")

pathways.show = c("FN1", "ncWNT", "EPHA", "NOTCH")

similarity_matrix_lift = data.frame(matrix(nrow = 4, ncol = length(pathways.show)+1))
colnames(similarity_matrix_lift) = c("Time",pathways.show)
similarity_matrix_lift$Time = c(6,9,12,15)

similarity_matrix_lift[1,2:ncol(similarity_matrix_lift)] = d6similarity[pathways.show,2]
similarity_matrix_lift[2,2:ncol(similarity_matrix_lift)] = d9similarity[pathways.show,2]
similarity_matrix_lift[3,2:ncol(similarity_matrix_lift)] = d12similarity[pathways.show,2]
similarity_matrix_lift[4,2:ncol(similarity_matrix_lift)] = d15similarity[pathways.show,2]

norm_sim_mat_lift = similarity_matrix_lift/max(similarity_matrix_lift[,2:ncol(similarity_matrix_lift)]) #normalize to max value of signaling pathways
norm_sim_mat_lift$Time = similarity_matrix_lift$Time #revert time to proper days again

#save(similarity_matrix_lift, norm_sim_mat_lift, file = "/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/similaritymatrices_LIFTUP.RData")

# # Specify id.vars: the variables to keep but not split apart on
# df = melt(norm_sim_mat, id.vars=c("Time"))

ggplot(data=melt(norm_sim_mat, id.vars=c("Time")), aes(x=Time, y=value, fill=variable)) +
  geom_bar(stat="identity",  position=position_dodge()) + theme_classic()

#Cell Type PCP Signaling Per Day

#load in individual cellchat objects that have been processed individually
D6ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D6ALLDB_lacZ.RDS")
D9ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D9ALLDB_lacZ.RDS")
D12ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D12ALLDB_lacZ.RDS")
D15ALLDB_lacZ <- readRDS("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/CellChat/D15ALLDB_lacZ.RDS")

sig = "NOTCH"
pat = "incoming"
netAnalysis_signalingRole_heatmap_topbar(D6ALLDB_lacZ, pattern = pat, signaling = sig, width = 5, height = 6)
netAnalysis_signalingRole_heatmap_topbar(D9ALLDB_lacZ, pattern = pat, signaling = sig, width = 5, height = 6)
netAnalysis_signalingRole_heatmap_topbar(D12ALLDB_lacZ, pattern = pat, signaling = sig, width = 5, height = 6)
netAnalysis_signalingRole_heatmap_topbar(D15ALLDB_lacZ, pattern = pat, signaling = sig, width = 5, height = 6)

#read in tables of incoming/outgoing signaling values per day
incomingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/ncWNT_incomingstrength_perday.csv")
outgoingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/ncWNT_outgoingstrength_perday.csv")

incomingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/FN1_incomingstrength_perday.csv")
outgoingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/FN1_outgoingstrength_perday.csv")

incomingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/EPHA_incomingstrength_perday.csv")
outgoingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/EPHA_outgoingstrength_perday.csv")

incomingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/NOTCH_incomingstrength_perday.csv")
outgoingdays <- read_csv("AQPhan/scRNA-seq/iPSC_Reprog/CellChat/NOTCH_outgoingstrength_perday.csv")

incomingdays$Type = "Incoming Signaling"
outgoingdays$Type = "Outgoing Signaling"

signaling = melt(rbind(incomingdays,outgoingdays), id.vars = c("Day", "Type"))

ggplot(signaling[signaling$Day==6,], aes(x=variable, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity", position = position_dodge()) + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45))
ggplot(signaling[signaling$Day==9,], aes(x=variable, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity", position = position_dodge()) + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45))
ggplot(signaling[signaling$Day==12,], aes(x=variable, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity", position = position_dodge()) + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45))
ggplot(signaling[signaling$Day==15,], aes(x=variable, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity", position = position_dodge()) + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45))


ggplot(signaling[signaling$Type=="Incoming Signaling",], aes(x=Day, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity") + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
ggplot(signaling[signaling$Type=="Outgoing Signaling",], aes(x=Day, y=value, group=Type, fill=variable, alpha = Type)) + geom_bar(stat = "identity") + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ggplot(signaling, aes(x=Day, y=value, group=Day, fill=variable, alpha = Type)) + geom_bar(stat = "identity", position = position_dodge()) + scale_alpha_manual(values=c(1, 0.5)) + theme_classic() + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
