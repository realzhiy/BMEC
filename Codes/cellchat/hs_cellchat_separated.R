library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(multtest)
library(metap)

setwd("~/Desktop/Applied_Genomics/BMEC/symphony")

# load human seurat object
human <- readRDS("human_combined_annotated1.rds")

t2dm <- subset(human, subset = orig.ident == "PDC05")

# Subset Seurat object for "cardiac" containing PDC34, PDC87, and PDC22
cardiac <- subset(human, subset = orig.ident %in% c("PDC34", "PDC87", "PDC22"))

# Subset Seurat object for "noncardiac" containing PDC91, C36, and C48
noncardiac <- subset(human, subset = orig.ident %in% c("PDC91", "C36", "C48"))



# Cell chat
## Create cellchat object
DefaultAssay(t2dm) <- "RNA"
cc.t2d <- createCellChat(t2dm, group.by = "ident")

# Setting ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB# simply use the default CellChatDB
# set the used database in the object
cc.t2d@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cc.t2d <- subsetData(cc.t2d) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cc.t2d <- identifyOverExpressedGenes(cc.t2d)
cc.t2d <- identifyOverExpressedInteractions(cc.t2d)

## Compute the communication probability and infer cellular communication network
cc.t2d <- computeCommunProb(cc.t2d)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cc.t2d <- filterCommunication(cc.t2d, min.cells = 10)

# a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cc.t2d)

# NB: The inferred intercellular communication network of each ligand-receptor 
# pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cc.t2d <- computeCommunProbPathway(cc.t2d)

## Calculate the aggregated cell-cell communication network
cc.t2d <- aggregateNet(cc.t2d)

cc.t2d <- netAnalysis_computeCentrality(cc.t2d)

pdf("t2d_circle_interaction_countandweight.pdf")
groupSize <- as.numeric(table(cc.t2d@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.t2d@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.t2d@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("t2d_edgeweights.pdf")
mat <- cc.t2d@net$weight
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2), xpd = TRUE)  # Adjusted plot margins
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


saveRDS(cc.t2d, file = "cellchat_t2d.rds")



## Create cellchat object
DefaultAssay(cardiac) <- "RNA"
cc.cardiac <- createCellChat(cardiac, group.by = "ident")

# Setting ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB# simply use the default CellChatDB
# set the used database in the object
cc.cardiac@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cc.cardiac <- subsetData(cc.cardiac) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cc.cardiac <- identifyOverExpressedGenes(cc.cardiac)
cc.cardiac <- identifyOverExpressedInteractions(cc.cardiac)

## Compute the communication probability and infer cellular communication network
cc.cardiac <- computeCommunProb(cc.cardiac)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cc.cardiac <- filterCommunication(cc.cardiac, min.cells = 10)

# a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cc.cardiac)

# NB: The inferred intercellular communication network of each ligand-receptor 
# pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cc.cardiac <- computeCommunProbPathway(cc.cardiac)

## Calculate the aggregated cell-cell communication network
cc.cardiac <- aggregateNet(cc.cardiac)

cc.cardiac <- netAnalysis_computeCentrality(cc.cardiac)

pdf("cardiac_circle_interaction_countandweight.pdf")
groupSize <- as.numeric(table(cc.cardiac@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.cardiac@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.cardiac@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("cardiac_edgeweights.pdf")
mat <- cc.cardiac@net$weight
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2), xpd = TRUE)  # Adjusted plot margins
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

saveRDS(cc.cardiac, file = "cellchat_cardiac.rds")



## Create cellchat object
DefaultAssay(noncardiac) <- "RNA"
cc.noncardiac <- createCellChat(noncardiac, group.by = "ident")

# Setting ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB# simply use the default CellChatDB
# set the used database in the object
cc.noncardiac@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cc.noncardiac <- subsetData(cc.noncardiac) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cc.noncardiac <- identifyOverExpressedGenes(cc.noncardiac)
cc.noncardiac <- identifyOverExpressedInteractions(cc.noncardiac)

## Compute the communication probability and infer cellular communication network
cc.noncardiac <- computeCommunProb(cc.noncardiac)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cc.noncardiac <- filterCommunication(cc.noncardiac, min.cells = 10)

# a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cc.noncardiac)

# NB: The inferred intercellular communication network of each ligand-receptor 
# pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cc.noncardiac <- computeCommunProbPathway(cc.noncardiac)

## Calculate the aggregated cell-cell communication network
cc.noncardiac <- aggregateNet(cc.noncardiac)

cc.noncardiac <- netAnalysis_computeCentrality(cc.noncardiac)

pdf("noncardiac_circle_interaction_countandweight.pdf")
groupSize <- as.numeric(table(cc.noncardiac@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.noncardiac@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.noncardiac@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("noncardiac_edgeweights.pdf")
mat <- cc.noncardiac@net$weight
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2), xpd = TRUE)  # Adjusted plot margins
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

saveRDS(cc.noncardiac, file = "cellchat_noncardiac.rds")



