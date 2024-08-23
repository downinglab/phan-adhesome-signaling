# Supplemental Figure 7 Monocle3 Pseudotime Trajectory Analysis
# 10/5/21
# By Andrew Phan

library(Seurat)
library(monocle3)
library(SeuratWrappers)

rm(list=ls()) #clears variables
cat("\014") #clears console

#combined <- readRDS("~/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/A549_combined_reclustered.rds") #load in A549 Seurat data
load(file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/Monocle3/Saved_objs/completed_A549_monocle3+PCP.RDATA") #load in finished analysis with 9 clusters

#plot a single gene signature score over pseudotime with a loess
#Input: dataframe with pseudotime and gene signature values; the column name for the x axis; the column name for the y axis; the number of breaks for the bins on the x axis [NOTE: TOO MANY BREAKS WILL PREVENT LOESS LINE FROM PLOTTING CORRECTLY]; a loess span value [Default=0.2]
#Output: A plot
plot_GeneSigScore <- function(df, x = "Pseudotime", y, breaks=50, span=0.2, ylim = NA){
  dfcut = df
  dfcut[[x]] <- cut(dfcut[[x]], breaks = breaks) #find breaks
  dfagg=aggregate(dfcut, by=list(dfcut[[x]]), mean) #average value per break

  lw1 <- loess(df[[y]] ~ df[[x]], data=df, span = span)

  if(!is.na(ylim)) #if ylim given, plot using limits
    plot(dfagg[[y]] ~ Group.1, data=dfagg,pch=19,cex=0.1, ylim = ylim)
  else
    plot(dfagg[[y]] ~ Group.1, data=dfagg,pch=19,cex=0.1)
  # pl <- ggplot(data = dfagg, aes(x=dfagg[[y]],y=dfagg[[y]])) + geom_point()
  # print(pl)
  j <- order(df[[x]])
  lines(df[[x]][j],lw1$fitted[j],col="red",lwd=3)
}

DimPlot(combined, reduction = "umap")
#separate conditions
CSRP1seurat = subset(x = combined, subset = orig.ident == "A549_D8_CSRP1_KD")
lacZseurat = subset(x = combined, subset = orig.ident == "A549_D8_LacZ")

seurat_obj = lacZseurat
seurat_obj = CSRP1seurat

cds <- as.cell_data_set(seurat_obj)

cds <- cluster_cells(cds)
#Retrieve the clustering from the Seurat object
#cds@clusters$UMAP$clusters <- seurat_obj@meta.data$seurat_clusters
#names(cds@clusters$UMAP$clusters) <- Cells(seurat_obj) #rename the cells
cds@clusters$UMAP$clusters <- seurat_obj@active.ident

plot_cells(cds, color_cells_by = "cluster")

#cds <- learn_graph(cds, learn_graph_control = list(ncenter = 550))
cds <- learn_graph(cds)

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")

cds <- order_cells(cds)

plot_cells(cds, #see results of pseudotime
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

############ a helper function to identify the root principal points: ##############
get_earliest_principal_node <- function(cds){ #finds vertex with most cells belonging to it given condition
  cell_ids <- which(colData(cds)[, "orig.ident"] == "A549_D8_CSRP1_KD")

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#################

traj_lacZ_TROPH <- choose_graph_segments(cds) #subset cells from a trajectory

cds_sub = traj_s3_NCC

cds_sub <- preprocess_cds(cds_sub, num_dim = 50)

cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP")

cds_sub <- cluster_cells(cds_sub)
cds_sub@clusters$UMAP$clusters <- CSRP1seurat@active.ident

plot_cells(cds_sub, group_cells_by = "cluster")

cds_sub <- learn_graph(cds_sub)

plot_cells(cds_sub,
           color_cells_by = "Time",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds_sub <- order_cells(cds_sub)

plot_cells(cds_sub, color_cells_by = "pseudotime")


cds_sub2 = cds_sub
cds_sub2@clusters$UMAP$clusters <- CSRP1seurat@active.ident

plot_cells(cds_sub2, group_cells_by = "cluster")


cds_subset <- choose_cells(cds)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=64)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < .1))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.1)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module,
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

###############

########### plot gene signatures over pseudotime ##############

#extract pseudotime values and signature scores

CSRP1seurat <- AddMetaData(
  object = CSRP1seurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

lacZseurat <- AddMetaData(
  object = lacZseurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

#need to update by running this line
seurat_obj = lacZseurat

#df = as.data.frame(seurat_obj@meta.data$orig.ident)
#colnames(df) = "Condition"
#df$Pseudotime = seurat_obj@meta.data$Pseudotime

#reorder and extract pseudotime from df object
df = df_lacZ[order(as.numeric(row.names(df_lacZ))), 1:2]

df$clust = seurat_obj@meta.data$Clusters
df$CSC.SIG = seurat_obj@meta.data$CSC1
df$EMT = seurat_obj@meta.data$EMT.LUAD1
df$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY_A5491


df <- df[order(df$Pseudotime),] #order based on pseudotime

seurat_obj = CSRP1seurat

#df2 = as.data.frame(seurat_obj@meta.data$orig.ident)
#colnames(df2) = "Condition"
#df2$Pseudotime = seurat_obj@meta.data$Pseudotime

#reorder and extract pseudotime from df object
df2 = df_CSRP1[ order(as.numeric(row.names(df_CSRP1))), 1:2]

df2$clust = seurat_obj@meta.data$Clusters
df2$CSC.SIG = seurat_obj@meta.data$CSC1
df2$EMT = seurat_obj@meta.data$EMT.LUAD1
df2$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY_A5491

df2 <- df2[order(df2$Pseudotime),] #order based on pseudotime

matplot(x=df$Pseudotime, y=df[,3:ncol(df)], type = c("b")) #matplot can plot multiple columns in one graph and add lines

########### find z-score = (value-mean)/SD ##############
#did not use, optional
dfzscore = df
for(x in 3:ncol(dfzscore)){ #for every signature score (column), calculate the z score
  mean = mean(df[,x])
  sd = sd(df[,x])
  dfzscore[,x] = (dfzscore[,x]-mean)/sd
}

dfzscore = dfzscore[!is.infinite(dfzscore$Pseudotime),]

plot_GeneSigScore(df, y="PCP", breaks = 15, span = 0.4) #CSRP1 CSC.sig zscore
##########

df = df[!is.infinite(df$Pseudotime),] #remove infinite pseudotime cells
df2 = df2[!is.infinite(df2$Pseudotime),] #remove infinite pseudotime cells

sig = "PCP.SIG"
plot_GeneSigScore(df, y=sig, breaks = 10, span = 0.2, ylim = c(-0.02, 0.025))
plot_GeneSigScore(df2, y=sig, breaks = 10, span = 0.2, ylim = c(-0.02, 0.025))


#FIGURES:
par(mfrow=c(1,2))    # set the plotting area into a 1x2 array
plot_GeneSigScore(df_lacZ, y=sig, breaks = 15, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
plot_GeneSigScore(df_CSRP1, y=sig, breaks = 15, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
#save as 10x20in pdf

##################################

FeaturePlot(lacZseurat, c("CSC1", "EMT.LUAD1"), pt.size = 0.1) & scale_color_viridis_c()

x="Pseudotime"
breaks = 20

dfcut = df
dfcut[[x]] <- cut(dfcut[[x]], breaks = breaks) #find breaks
dfagg=aggregate(dfcut, by=list(dfcut[[x]]), mean) #average value per break
plot(dfagg[["PCP.SIG"]] ~ Group.1, data=dfagg,pch=19,cex=0.1)

lw1 <- loess(formula = PCP.SIG ~ Pseudotime, data = df)
plot(PCP.SIG ~ Pseudotime, data=df,pch=19,cex=0.1)
j <- order(df$Pseudotime)
lines(df$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)

lw1 <- loess(formula = PCP.SIG ~ Pseudotime, data = df2)
plot(PCP.SIG ~ Pseudotime, data=df2,pch=19,cex=0.1)
j <- order(df2$Pseudotime)
lines(df2$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)

########################## CellChat
rm(list=ls()) #clears variables
cat("\014") #clears console

library(CellChat)
library(patchwork)
library(future)
library(gplots)
library(viridis)

########### Fixed Compute Similarity Related Functions ############
# source - https://github.com/sqjin/CellChat/issues/167
# Copied + pasted with no due diligence.

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
  Y <- runUMAP(Similarity, min.dist = 0.3, n.neighbors = k)
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

#load in Seurat object of combined lacZ and CSRP1 KD scRNA-seq
combined <- readRDS("~/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/A549_combined_reclustered.rds")
combined <- RenameIdents(object = combined, `0` = "1", `1` = "2", `2` = "3", `3` = "4") #rename because can't include `0`

#LacZcellchat <- readRDS("~/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/")
#CSRP1cellchat <- readRDS("~/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/")
# object.list <- list(L = LacZcellchat, S = CSRP1cellchat)
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# cellchat

############ INDIVIDUAL CELLCHAT ANALYSIS ####################
########### Convert Seurat object to CellChat object ############

library(Seurat) #need Seurat package loaded for SplitObject and GetAssayData

splitSeurat = SplitObject(combined, split.by = "orig.ident")

#create LacZ CellChat object
data.input <- GetAssayData(object = splitSeurat$A549_D8_LacZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(splitSeurat$A549_D8_LacZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#LacZcellchat <- createCellChat(object = splitSeurat$LacZ, assay = "RNA", meta = meta, group.by = "group")
LacZcellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#create SHROOM3 CellChat object
data.input <- GetAssayData(object = splitSeurat$A549_D8_CSRP1_KD, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(splitSeurat$A549_D8_CSRP1_KD)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

CSRP1cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")


######### Set the ligand-receptor interaction database ############

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

#CellChatDB.use <- subsetDB(CellChatDB, search = c("ECM-Receptor", "Cell-Cell Contact")) #focus on cell-ECM and cell-cell interactions (not para/autocrine)
CellChatDB.use <- subsetDB(CellChatDB, search = c("ECM-Receptor", "Cell-Cell Contact", "Secreted Signaling")) #focus on cell-ECM and cell-cell interactions (not para/autocrine)
#CellChatDB.use = CellChatDB #use all cell signaling pathways
LacZcellchat@DB <- CellChatDB.use
CSRP1cellchat@DB <- CellChatDB.use

############ Preprocessing the expression data for cell-cell communication analysis ############

LacZcellchat <- subsetData(LacZcellchat) # This step is necessary even if using the whole database
CSRP1cellchat <- subsetData(CSRP1cellchat) # This step is necessary even if using the whole database

future::plan("multiprocess", workers = 4) # do parallel

LacZcellchat <- identifyOverExpressedGenes(LacZcellchat)
LacZcellchat <- identifyOverExpressedInteractions(LacZcellchat)

CSRP1cellchat <- identifyOverExpressedGenes(CSRP1cellchat)
CSRP1cellchat <- identifyOverExpressedInteractions(CSRP1cellchat)

#cellchat <- projectData(cellchat, PPI.human) # project gene expression data onto PPI network (optional)

############ Inference of cell-cell communication network ##########

LacZcellchat <- computeCommunProb(LacZcellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
LacZcellchat <- filterCommunication(LacZcellchat, min.cells = 10)

CSRP1cellchat <- computeCommunProb(CSRP1cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
CSRP1cellchat <- filterCommunication(CSRP1cellchat, min.cells = 10)

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

#Infer the cell-cell communication at a signaling pathway level
LacZcellchat <- computeCommunProbPathway(LacZcellchat)
CSRP1cellchat <- computeCommunProbPathway(CSRP1cellchat)

#Calculate the aggregated cell-cell communication network
LacZcellchat <- aggregateNet(LacZcellchat)
CSRP1cellchat <- aggregateNet(CSRP1cellchat)

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
CSRP1cellchat@netP$pathways

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
CSRP1cellchat <- netAnalysis_computeCentrality(CSRP1cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
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

# #Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "incoming")
ht1 + ht2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(CSRP1cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(CSRP1cellchat, pattern = "incoming")
ht1 + ht2

# #Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# ht1 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "outgoing")
# ht2 <- netAnalysis_signalingRole_heatmap(LacZcellchat, pattern = "incoming")
# ht1 + ht2

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)
library(ggalluvial)

#NOTE: NEED TO ADJUST nPatterns PER CONDITION
selectK(LacZcellchat, pattern = "outgoing")
nPatterns = 3 #change to before where both graphs drop off together
LacZcellchat <- identifyCommunicationPatterns(LacZcellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(LacZcellchat, pattern = "outgoing")

selectK(CSRP1cellchat, pattern = "outgoing")
nPatterns = 3
CSRP1cellchat <- identifyCommunicationPatterns(CSRP1cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(CSRP1cellchat, pattern = "outgoing")

selectK(LacZcellchat, pattern = "incoming")
nPatterns = 2
LacZcellchat <- identifyCommunicationPatterns(LacZcellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(LacZcellchat, pattern = "incoming")

selectK(CSRP1cellchat, pattern = "incoming")
nPatterns = 2
CSRP1cellchat <- identifyCommunicationPatterns(CSRP1cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(CSRP1cellchat, pattern = "incoming")

#Manifold and classification learning analysis of signaling networks
LacZcellchat <- computeNetSimilarity(LacZcellchat, type = "functional")
LacZcellchat <- netEmbedding(LacZcellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
LacZcellchat <- netClustering(LacZcellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(LacZcellchat, type = "functional", label.size = 3.5)

LacZcellchat <- computeNetSimilarity(LacZcellchat, type = "structural")
LacZcellchat <- netEmbedding(LacZcellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
LacZcellchat <- netClustering(LacZcellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(LacZcellchat, type = "structural", label.size = 3.5)

CSRP1cellchat <- computeNetSimilarity(CSRP1cellchat, type = "functional")
CSRP1cellchat <- netEmbedding(CSRP1cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
CSRP1cellchat <- netClustering(CSRP1cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(CSRP1cellchat, type = "functional", label.size = 3.5)

CSRP1cellchat <- computeNetSimilarity(CSRP1cellchat, type = "structural")
CSRP1cellchat <- netEmbedding(na.omit(CSRP1cellchat), type = "structural") #add na.omit() if giving error: Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)
#> Manifold learning of the signaling networks for a single dataset
CSRP1cellchat <- netClustering(CSRP1cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(CSRP1cellchat, type = "structural", label.size = 3.5)

#saveRDS(LacZcellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/LacZcellchat_individual_A549.rds")
#saveRDS(CSRP1cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChatCSRP1cellchat_individual_A549.rds")

##################### Merge the LacZ and S3 CellChat Objects After Individual Analysis #######################
#create a merged CellChat object of LacZ and S3 conditions
object.list <- list(lacZ = LacZcellchat, CSRP1 = CSRP1cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

save(cellchat, LacZcellchat, CSRP1cellchat, file="/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChatLacZ+S3cellchat_combined+allsignaling.RDATA")
#saveRDS(S3cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/CellChat/LacZ+S3cellchat_combined.rds")

############### COMPARATIVE CELLCHAT ANALYSIS ###############
############### Predict general principles of cell-cell communication ###############
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#big difference in interaction strength

#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#red/blue colored edges represent increased/decreased signaling in the 2nd dataset compared to the 1st

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
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
#group 1 is much more important in CSRP1 KD for incoming and outgoing interactions

############## Identify the conserved and context-specific signaling pathways ###############
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

#Error in data.frame(x = Y[, 1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),  :
#                      arguments imply differing number of rows: 60, 0

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
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
heatmap.2(t(rank), Colv = NA, Rowv = NA, col = magma(256)) #heatmap of log(distance) of pathways LvC

#Identify and visualize the conserved and context-specific signaling pathways
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)

# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Outgoing: Similar FN1 and NOTCH, less EPHA

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#incoming: higher FN1, less NOTCH and EPHA

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Overall: similar FN1, less EPHA and NOTCH

#desmosome less in this and in iPSC reprogramming in KD conditions compared to controls

################### Identify the upregulated and down-regulated signaling ligand-receptor pairs ######################
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:4),  comparison = c(1, 2), angle.x = 45)
#FN1 signaling mediated by CD44 (CSC marker, allows for metastasis/EMT). Makes sense in this context

gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CSRP1", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in CSRP1", angle.x = 45, remove.isolate = T)
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
net.down <- subsetCommunication(cellchat, net = net, datasets = "CSRP1",ligand.logFC = -0.1, receptor.logFC = -0.1)

#see individual signaling genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:11), targets.use = c(1:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(1:4), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(1:4), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#Error in netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11),  :
#   No signaling links are inferred!

#saveRDS(cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/combinedCellchat_A549_processed.rds")

############### Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram ####################
#edge diagram
pathways.show <- c("NOTCH")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#heatmap
pathways.show <- c("EPHA")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
#no difference in EPHA in cluster 2 (CSCs)
#less NOTCH incoming to cluster 2

# Chord diagram
pathways.show <- c("DESMOSOME")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

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

pathways.show <- c("FN1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(1:4), targets.use = c(1:4), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("EPHA")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(1:4), targets.use = c(1:4), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("NOTCH")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, sources.use = c(1:4), targets.use = c(1:4), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], signaling = "EPHA", sources.use = c(1), targets.use = c(1), slot.name = "netP", title.name = paste0("Sig", names(object.list)[i]), legend.pos.x = 10)
}

#Vertex FIGURES:
#visualize vertex graph with clusters 1 and 2 on left
i=1 #choose which condition (LacZ)
i=2 #CSRP1
vertex.receiver = c(1:2)
pathways.show = "NOTCH"
netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
#Notch decreases in CSRP1

i=1 #choose which condition (LacZ)
i=2 #CSRP1
vertex.receiver = c(1:2)
pathways.show = "FN1"
netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
#FN1 connectivity seems to be the same, but bigger sources for 1+2

i=1 #choose which condition (LacZ)
i=2 #CSRP1
vertex.receiver = c(1:2)
pathways.show = "EPHA"
netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))


# Aggregated Chord diagram
#LacZ in
# 1: 9
# 2: 2, 3, 5, 6
# 3: 1, 4, 8, 10
group.cellType <- c(1,1,2,3) # grouping cell clusters based on LacZ in
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("FN1")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, remove.isolate = T, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

################# Compare the signaling gene expression distribution between different datasets ################
# cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("L", "S")) # set factor level
plotGeneExpression(cellchat, signaling = "FN1", split.by = "orig.ident", colors.ggplot = T, split.plot = T)

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
vertex.receiver = c(1,2)
pathways.show = "FN1"
netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))

#save(LacZcellchat, CSRP1cellchat, cellchat, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/A549/CellChat/indiv+combined_processed_cellchat.RDATA")
