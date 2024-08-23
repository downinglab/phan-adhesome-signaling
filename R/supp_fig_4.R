# Supplemental Figure 4 Chemically Induced Pluripotent Cells CellChat
# By Nolan Origer

# Libraries and constants ---------------------------------------------------------------
# https://www.nature.com/articles/s41586-022-04593-5
# hCiPs = human Chemically induced Pluripotent Cells

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(DoMultiBarHeatmap)
library(ggpubr)

stage_cols <- c(hADSCs = "#1a7ad6", S1D0 = "#84c4f8", 
                S1D0.5 = "#cbe7cc", S1D1 = "#8bcb8d", S1D2 = "#58b55b", S1D4 = "#409a44", S1D8 = "#2c7a30",
                S2D4 = "#feddab", S2D8 = "#feb03d", S2D12 = "#fc9000", S2D20 = "#ef6e00",
                S3D4 = "#f7b3b5", S3D8 = "#f34940", S3D12 = "#cc2b2b",
                S4D1 = "#deb8e5", S4D2 = "#b35ac3", S4D4 = "#9323ac", S4D10 = "#6c1b9b",
                hCiPSCs = "#3241a1", H1 = "#959fd6")

traj_cols <- c(R0 = "#1a7ad6", R1 = "#84c4f8", 
               R2 = "#cbe7cc", R3 = "#8bcb8d", R4 = "#58b55b", R5 = "#409a44", R6 = "#2c7a30",
               R7 = "#feddab", R8 = "#fc9000", R9 = "#ef6e00",
               R10 = "#f7b3b5", R11 = "#cc2b2b",
               R12 = "#deb8e5", R13 = "#9323ac",
               Failed = "#e5e5e5",
               hCiPSCs = "#3241a1", H1 = "#959fd6")
# hCiP Data ----------------------------------------------------------------

### Define hCiP sample names
cip_times <- c("hADSCs",
               "S1D0", "S1D0.5", "S1D1", "S1D2", "S1D4", "S1D8", "S1D16",
               "S2D4", "S2D8", "S2D12", "S2D20", "S2D24", 
               "S3D4", "S3D8", "S3D12", 
               "S4D1", "S4D2", "S4D4", "S4D10",
               "hCiPSCs", "H1")

### Downloaded files from GSE178325 and renamed folders according to scheme above
### Now load in Seurat objects for each sample/timepoint
for (i in cip_times){
  cat(paste0("Processing sample ",i,"\n"))
  assign(paste0("cip_data_", i), Read10X(data.dir = paste0("###DIRECTORY###/", i)))
  assign(paste0("cip_", i), CreateSeuratObject(counts = eval(parse(text = paste0("cip_data_", i))), project = i))
}

### Ensure no cell name overlap
for (i in 1:length(cip_times)){
  cat(paste0("Processing object ", i,"\n"))
  assign(paste0("cip_", cip_times[i]), RenameCells(object = eval(parse(text = paste0("cip_", cip_times[i]))), add.cell.id = cip_times[i]))
}

### Combine the objects
cip <- merge(cip_hADSCs, y = c(cip_S1D0, cip_S1D0.5, cip_S1D1, cip_S1D2, cip_S1D4, cip_S1D8, cip_S1D16,
                               cip_S2D4, cip_S2D8, cip_S2D12, cip_S2D20, cip_S2D24,
                               cip_S3D4, cip_S3D8, cip_S3D12,
                               cip_S4D1, cip_S4D2, cip_S4D4, cip_S4D10,
                               cip_hCiPSCs, cip_H1))

### Store stage name as new metadata column
Idents(cip) <- "orig.ident"
cip$stage <- cip$orig.ident

### Save this section
saveRDS(cip, "hcip_raw.RDS")

# hCiP basic quality control -----------------------------------------

### Load previous section
cip <- readRDS("hcip_raw.RDS")

### Load annotations from source paper github
annotations <- read.delim("Fig2b_anno.tsv", header = T)

### Keep cells based on paper annotations
cip <- subset(cip, cells = rownames(annotations))

### Add metadata from paper annotations
cip <- AddMetaData(cip, metadata = annotations)
rm(annotations)

### Factor stages
cip$stage <- factor(cip$stage, levels = c("hADSCs",
                                          "S1D0", "S1D0.5", "S1D1", "S1D2", "S1D4", "S1D8", "S1D16",
                                          "S2D4", "S2D8", "S2D12", "S2D20", "S2D24", 
                                          "S3D4", "S3D8", "S3D12", 
                                          "S4D1", "S4D2", "S4D4", "S4D10",
                                          "hCiPSCs", "H1"))

### Process as done in paper according to source github
DefaultAssay(cip) <- "RNA"
cip <- JoinLayers(cip)
cip <- NormalizeData(cip)
cip <- FindVariableFeatures(cip)
cip <- ScaleData(cip)
cip <- RunPCA(cip)
cip <- FindNeighbors(cip, dims = 1:20)
cip <- FindClusters(cip, resolution = 1)
cip <- RunUMAP(cip, dims = 1:20)

### Plot FA according to coordinates from source paper
ggplot(cip@meta.data,aes(x=fa_1,y=fa_2,color=stage)) + 
  geom_point(size = 0.1) + 
  theme_classic() +
  labs(title = NULL, x = "FA1", y = "FA2") +
  guides(x = ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),
                                         trunc_upper = unit(2.9, "cm")), 
         y = ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),
                                         trunc_upper = unit(2.9, "cm"))) +
  guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.25, "cm"), 
                                               type = "closed")),
        axis.title = element_text(hjust = 0)) +
  scale_color_manual(name = NULL,
                     values = stage_cols) + 
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)

### Rename trajectory labels for readability
Idents(cip) <- "traj"
cip <- RenameIdents(cip,
                    'R_0' = "R0",
                    'R_1' = "R1",
                    'R_2' = "R2",
                    'R_3' = "R3",
                    'R_4' = "R4",
                    'R_5' = "R5",
                    'R_6' = "R6",
                    'R_7' = "R7",
                    'R_8' = "R8",
                    'R_9' = "R9",
                    'R_10' = "R10",
                    'R_11' = "R11",
                    'R_12' = "R12",
                    'R_13' = "R13",
                    'F' = "Failed",
                    'prime' = "hCiPSCs",
                    'ES' = "H1")
cip$traj <- Idents(cip)

### Plot FA colored by trajectory
ggplot(cip@meta.data,aes(x=fa_1,y=fa_2,color=traj)) + 
  geom_point(size = 0.1) + 
  theme_classic() +
  labs(title = NULL, x = "FA1", y = "FA2") +
  guides(x = ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),
                                         trunc_upper = unit(2.9, "cm")), 
         y = ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),
                                         trunc_upper = unit(2.9, "cm"))) +
  guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.25, "cm"), 
                                               type = "closed")),
        axis.title = element_text(hjust = 0)) +
  scale_color_manual(name = NULL,
                     values = traj_cols) + 
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)

### Save this section
saveRDS(cip, "hcip_qc.RDS")

# hCiP module scoring ------------------------------------------------

### Load previous section
cip <- readRDS("hcip_qc.RDS")

### Load CellChat database to get genes from each signaling pathway
library(CellChat)
### Get EPHA genes and process into gene list
score.EPHA <- subsetDB(CellChatDB.human, search = c("EPHA"), key = "pathway_name") # use Secreted Signaling
score.EPHA <- unique(c(score.EPHA$interaction$ligand, score.EPHA$interaction$receptor)) %>% 
  strsplit("_") %>%  unlist() %>% unique() %>% list()

### Get FN1 genes and process into gene list
score.FN1 <- subsetDB(CellChatDB.human, search = c("FN1"), key = "pathway_name") # use Secreted Signaling
score.FN1 <- unique(c(score.FN1$interaction$ligand, score.FN1$interaction$receptor)) %>% 
  strsplit("_") %>%  unlist() %>% unique() %>% list()

### Get NOTCH genes and process into gene list
score.NOTCH <- subsetDB(CellChatDB.human, search = c("NOTCH"), key = "pathway_name") # use Secreted Signaling
score.NOTCH <- unique(c(score.NOTCH$interaction$ligand, score.NOTCH$interaction$receptor)) %>% 
  strsplit("_") %>%  unlist() %>% unique() %>% list()

### Add module score using above gene lists
cip <- AddModuleScore(cip, features = score.EPHA, name = "EPHA_Score")
cip <- AddModuleScore(cip, features = score.FN1, name = "FN1_Score")
cip <- AddModuleScore(cip, features = score.NOTCH, name = "NOTCH_Score")

### Rename as necessary
cip$EPHA_Score <- cip$EPHA_Score1
cip$FN1_Score <- cip$FN1_Score1
cip$NOTCH_Score <- cip$NOTCH_Score1
cip$EPHA_Score1 <- NULL
cip$FN1_Score1 <- NULL
cip$NOTCH_Score1 <- NULL

### Save this section
saveRDS(cip, "hcip_scored.RDS")

# hCiP CellChat -----------------------------------------------------------
### Ensure necessary libraries are loaded
library(CellChat)
library(NMF)
library(ggalluvial)

### Load previous section
cip <- readRDS("hcip_scored.RDS")

### Define stages with more than one trajectory label
stages <- c("S1D4", "S1D8", "S2D4", "S2D8", "S2D12", "S2D20", 
            "S3D4", "S3D8", "S3D12", "S4D1", "S4D2", "S4D4", "S4D10")

### Generate and save CellChat object for each stage
for(i in stages){
  meta <- data.frame(stage = cip$stage, 
                     traj = cip$traj)
  
  meta <- subset(meta, subset = stage == i)
  
  cellchat <- subset(cip, subset = stage == i) %>% 
    GetAssayData(layer = "data") %>%
    createCellChat(meta = meta, group.by = "traj")
  
  cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
  
  ### Set the database in the object to only cell-cell and cell-ECM
  cellchat@DB <- subsetDB(CellChatDB.human, search = c("ECM-Receptor","Cell-Cell Contact"), key = "annotation")
  
  ### Subset data based on database of interest and identify DEGs
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, paste0(i, "_cellchat.RDS"))
}

### Load CellChat object for each stage, apply labels, and add to object.list
group.new <- c("Failed", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13")
group.new <- factor(group.new, levels = c("Failed", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13"))
object.list <- list()
for(i in stages){
  object.list[i] <- readRDS(paste0(i, "_cellchat.RDS"))
}

### Apply same labels to all objects in list
object.list[[1]] <- liftCellChat(object.list[[1]], group.new = c("R5", "Failed", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13"))
object.list[[2]] <- liftCellChat(object.list[[2]], group.new = c("R5", "R6", "Failed", "R7", "R8", "R9", "R10", "R11", "R12", "R13"))
object.list[[3]] <- liftCellChat(object.list[[3]], group.new = c("R5", "R6", "R7", "R8", "Failed", "R9", "R10", "R11", "R12", "R13"))
object.list[[4]] <- liftCellChat(object.list[[4]], group.new = c("R5", "R6", "R7", "R8", "Failed", "R9", "R10", "R11", "R12", "R13"))
object.list[[5]] <- liftCellChat(object.list[[5]], group.new = c("R5", "R6", "R7", "R8", "R9", "Failed", "R10", "R11", "R12", "R13"))
object.list[[6]] <- liftCellChat(object.list[[6]], group.new = c("R5", "R6", "R7", "R8", "R9", "Failed", "R10", "R11", "R12", "R13"))
object.list[[7]] <- liftCellChat(object.list[[7]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "Failed", "R13"))
object.list[[8]] <- liftCellChat(object.list[[8]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "Failed", "R13"))
object.list[[9]] <- liftCellChat(object.list[[9]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "Failed", "R13"))
object.list[[10]] <- liftCellChat(object.list[[10]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "Failed", "R13"))
object.list[[11]] <- liftCellChat(object.list[[11]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "Failed", "R13"))
object.list[[12]] <- liftCellChat(object.list[[12]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "Failed"))
object.list[[13]] <- liftCellChat(object.list[[13]], group.new = c("R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "Failed"))

### Apply analysis functions
object.list <- lapply(object.list, function (x) setIdent(x, ident.use = "traj", levels = group.new))
object.list <- lapply(object.list, function (x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function (x) aggregateNet(x))
object.list <- lapply(object.list, function (x) netAnalysis_computeCentrality(x, slot.name = "netP"))

### Merge objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

### Choose signaling pathways of interest based on Fig. 3c and Supplemental Fig. 3a,b
pathway.union <- c("NOTCH","FN1","EPHA","HSPG","DESMOSOME","EPHB","SEMA5")

### Custom function with code from netAnalysis_signalingRole_heatmap for getting matrix
signaling_heatmap_matrix <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP"){
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
  } else if (pattern == "incoming") {
    mat <- t(incoming)
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  
  return(mat)
  
}

### Get matrices for each day
mat.list <- list()
pathway.union <- c("NOTCH","FN1","EPHA","HSPG","DESMOSOME","EPHB","SEMA5")
for(i in names(object.list)){
  ### Get signaling matrix then remove columns that are all zeros
  mat.list[[i]] <- signaling_heatmap_matrix(object.list[[i]], signaling = pathway.union, pattern = "all") %>%
    data.frame() %>%
    select(where(~ sum(.) != 0))
  colnames(mat.list[[i]]) <- paste0(colnames(mat.list[[i]]), ".", i)
  mat.list[[i]] <- sweep(mat.list[[i]], 1L, apply(mat.list[[i]], 1, max), '/', check.margin = FALSE)
  mat.list[[i]] <- tibble::rownames_to_column(mat.list[[i]], var = "Row.names")
}

### Merge the generated matrix list
mat.merge <- purrr::reduce(mat.list, merge, by = c("Row.names")) %>%
  tibble::column_to_rownames(var = "Row.names")
mat.merge[is.na(mat.merge)] <- 0

### Save this section
saveRDS(mat.merge, "mat_merge.RDS")

# hCiP Pseudobulk ---------------------------------------------------------

### Load previous section
cip <- readRDS("hcip_scored.RDS")

### Aggregate SHROOM3 expression by trajectory and stage
cip.bulk <- AggregateExpression(cip, 
                                assays = "RNA", 
                                features = "SHROOM3",
                                group.by = c("traj", "stage"))$RNA %>%
  as.matrix()

### Make column names match cellchat matrices
colnames(cip.bulk) <- gsub("_", ".", colnames(cip.bulk))
rownames(cip.bulk) <- "SHROOM3"
cip.bulk <- cip.bulk[1, c("Failed.S1D4", "Failed.S1D8", 
                          "Failed.S2D4", "Failed.S2D8", "Failed.S2D12", "Failed.S2D20", 
                          "Failed.S3D4", "Failed.S3D8", "Failed.S3D12",  
                          "Failed.S4D1", "Failed.S4D2", "Failed.S4D4", 
                          "R5.S1D4", "R5.S1D8", 
                          "R6.S1D8", 
                          "R7.S2D4", "R7.S2D8", "R7.S2D12", "R7.S2D20",  
                          "R8.S2D4", "R8.S2D8", "R8.S2D20", 
                          "R9.S2D12", "R9.S2D20", 
                          "R10.S3D4", "R10.S3D8", "R10.S3D12", "R10.S4D1", "R10.S4D2", "R10.S4D4", "R10.S4D10",
                          "R11.S3D12", 
                          "R11.S4D1", "R11.S4D2", "R11.S4D4",  
                          "R12.S3D4", "R12.S3D8", "R12.S3D12", 
                          "R12.S4D1", "R12.S4D2", "R12.S4D10",
                          "R13.S4D4", "R13.S4D10")] %>% as.matrix() %>% t()

### Save this section
saveRDS(cip.bulk, "cip_bulk_S3.RDS")

# hCiP Figures ------------------------------------------------------------

##### Supplementary Fig. 4a
cip <- readRDS("hcip_scored.RDS")

Idents(cip) <- cip$stage

### Change "features" to desired score
VlnPlot(cip, features = "NOTCH_Score", 
        group.by = "stage", 
        idents = levels(cip$stage)[6:20], 
        split.by = "fate", 
        split.plot = T, 
        pt.size = 0) +
  ggpubr::stat_compare_means(label = "p.signif") +
  stat_summary(aes(color = subset(cip, stage %in% levels(cip$stage)[6:20])$fate), fun.y = median, geom='point', size = 10, shape = 95) +
  scale_fill_manual(name = "Reprogramming Trajectory",
                    values = c(Progressing = "darkolivegreen3",
                               Failed = "firebrick3")) +
  scale_color_manual(name = "Median",
                     values = c(Progressing = "darkolivegreen1",
                                Failed = "firebrick1"))

##### Supplementary Fig. 4b
cip.shroom3 <- readRDS("cip_bulk_S3.RDS")
mat.merge <- readRDS("mat_merge.RDS")

### Scale SHROOM3 data from 0 to 1
cip.shroom3 <- sweep(cip.shroom3, 1L, apply(cip.shroom3, 1, max), '/', check.margin = FALSE)

### Use curated list of trajectory/stage based on cell sample size
cip.shroom3 <- cip.shroom3[, c("Failed.S1D4", "Failed.S1D8", 
                               "Failed.S2D4", "Failed.S2D8", "Failed.S2D12", "Failed.S2D20", 
                               "Failed.S3D4", "Failed.S3D8", "Failed.S3D12",  
                               "Failed.S4D1", "Failed.S4D2", "Failed.S4D4", 
                               "R5.S1D4", 
                               "R6.S1D8", 
                               "R7.S2D4", 
                               "R8.S2D8", "R8.S2D20", 
                               "R9.S2D20", 
                               "R10.S3D4", "R10.S3D8", "R10.S3D12", "R10.S4D1", "R10.S4D2", "R10.S4D4",
                               "R11.S3D12", 
                               "R11.S4D1", "R11.S4D2", "R11.S4D4",  
                               "R12.S3D12", 
                               "R12.S4D1", "R12.S4D2", "R12.S4D10",
                               "R13.S4D10")]

mat.merge <- mat.merge[, c("Failed.S1D4", "Failed.S1D8", 
                           "Failed.S2D4", "Failed.S2D8", "Failed.S2D12", "Failed.S2D20", 
                           "Failed.S3D4", "Failed.S3D8", "Failed.S3D12",  
                           "Failed.S4D1", "Failed.S4D2", "Failed.S4D4", 
                           "R5.S1D4", 
                           "R6.S1D8", 
                           "R7.S2D4", 
                           "R8.S2D8", "R8.S2D20", 
                           "R9.S2D20", 
                           "R10.S3D4", "R10.S3D8", "R10.S3D12", "R10.S4D1", "R10.S4D2", "R10.S4D4",
                           "R11.S3D12", 
                           "R11.S4D1", "R11.S4D2", "R11.S4D4",  
                           "R12.S3D12", 
                           "R12.S4D1", "R12.S4D2", "R12.S4D10",
                           "R13.S4D10")]

Heatmap(mat.merge, 
        col = colorRampPalette(c("white", "#4fb589", "#03461e"))(5), 
        na_col = "white", 
        name = "Relative \nStrength", 
        top_annotation = columnAnnotation(Sum = anno_barplot(colSums(mat.merge)), 
                                          SHROOM3 = as.numeric(cip.shroom3), 
                                          col = list(SHROOM3 = colorRamp2(c(0, 0.5, 1), c("white", "#4fb589", "#03461e"))),
                                          annotation_name_side = "left",
                                          show_legend = c(FALSE)), 
        right_annotation = rowAnnotation(Sum = anno_barplot(rowSums(mat.merge))), 
        cluster_rows = T,
        show_row_dend = F, 
        cluster_columns = F, 
        row_names_side = "left", 
        row_names_rot = 0, 
        column_names_rot = 90)
