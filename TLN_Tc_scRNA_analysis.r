## ---------------------------
##
## scRNA analysis for A Tessellated Lymphoid Network Provides Whole-Body T cell Surveillance in Zebrafish
##
## Author: Yiran Hou
##
## Date Created: 2022-10-16
##

#***********************************
# Loading required libraries
#***********************************
library(Seurat) # main scRNA analysis package
library(SeuratDisk) # H5 format save/load
library(DoubletFinder)
library(ggplot2) # plotting
library(dittoSeq) # Seurat-compatible visualizations
library(scCustomize) # customized plotting

#***********************************
# Data loading and basic QC
#***********************************
#1. Create Seurat object
## Control
ctrlT.data <- Read10X("Control/control_filtered_feature_bc_matrix/")
ctrlT <- CreateSeuratObject(counts = ctrlT.data, project = "ctrl", min.cells = 4)
## Infected
infectT.data <- Read10X("Infected/infected_filtered_feature_bc_matrix/")
infeT <- CreateSeuratObject(counts = infectT.data, project = "infected", min.cells = 4)

#2. Selecting cells for further analysis
# Calculate mitochondria contribution
ctrlT[["percent.mt"]] <- PercentageFeatureSet(ctrlT, pattern = "^mt-")
infeT[["percent.mt"]] <- PercentageFeatureSet(infeT, pattern = "^mt-")
# Calculate ribosomal contribution 
ctrlT[["percent.ribo"]] <- PercentageFeatureSet(ctrlT, pattern = "^rp[sl]")
infeT[["percent.ribo"]] <- PercentageFeatureSet(infeT, pattern = "^rp[sl]")
# Proportion of hemoglobin genes
hemo.genes <- c("hbaa1","hbaa2","hbae1.1","hbae1.3","hbae3","hbae5","hbba1","hbba2","hbbe1.1","hbbe2","hbbe3","hbbe1.2","hbbe1.3")
total_cts_per_cell_ctrl <- colSums(ctrlT@assays$RNA@counts)
hemo.genes.ctrl <- match(hemo.genes, rownames(ctrlT))
hemo.genes.ctrl <- hemo.genes.ctrl[!is.na(hemo.genes.ctrl)]
ctrlT[["percent.hb"]] <- colSums(ctrlT@assays$RNA@counts[rownames(ctrlT)[hemo.genes.ctrl], ])/total_cts_per_cell_ctrl
total_cts_per_cell_infect <- colSums(infeT@assays$RNA@counts)
hemo.genes.infe <- match(hemo.genes, rownames(infeT))
hemo.genes.infe <- hemo.genes.infe[!is.na(hemo.genes.infe)]
infeT[["percent.hb"]] <- colSums(infeT@assays$RNA@counts[rownames(infeT)[hemo.genes.infe], ])/total_cts_per_cell_infect
# Subset criteria for both samples:
# nFeature > 200, mt < 5. Leave higher end filtering to doublet finding step.
ctrlT <- subset(ctrlT, subset = nFeature_RNA > 200 & percent.mt < 5)
infeT <- subset(infeT, subset = nFeature_RNA > 200 & percent.mt < 5)

#3. Cell cycle scoring
# Load cell cycle genes
fishCCgenes <- readLines(con = "regev_lab_cell_cycle_genes_asFish.txt") 
s.genes <- fishCCgenes[1:42]
g2m.genes <- fishCCgenes[43:96]
# Normalize data
ctrlT <- NormalizeData(ctrlT)
infeT <- NormalizeData(infeT)
# Cell cycle scoring
ctrlT <- CellCycleScoring(object = ctrlT, g2m.features = g2m.genes,
    s.features = s.genes)
infeT <- CellCycleScoring(object = infeT, g2m.features = g2m.genes,
    s.features = s.genes)

#4. Initial linear/non-linear dimensional reduction
# SCTransform in place of normalization, finding variable features and scaling data.
ctrlT = SCTransform(ctrlT, vars.to.regress = c("nFeature_RNA","percent.mt"))
infeT = SCTransform(infeT, vars.to.regress = c("nFeature_RNA","percent.mt"))
# Linear and non-linear dimensional reduction for control
ctrlT = RunPCA(ctrlT, npcs = 50)
ctrlT = RunUMAP(ctrlT, dims = 1:30)
# Linear and non-linear dimensional reduction for infected
infeT = RunPCA(infeT, npcs = 50)
infeT = RunUMAP(infeT, dims = 1:30)

#5. Doublet finding
# [Control] Use pN 0.25, pK 0.09 for doublet finding
nExp <- round(ncol(ctrlT) * 0.03) # estimate number of doublets. 3% multiplet rate for 3000-4000 cells recovered
ctrlT <- doubletFinder_v3(ctrlT, sct = TRUE, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30)
# [Control] Doublet removal. Subset name comes from doubletFinder
ctrlT.sub <- subset(ctrlT, subset = DF.classifications_0.25_0.09_97 =="Singlet")
# [Infected] Use pN 0.25, pK 0.09 for doublet finding
nExp <- round(ncol(infeT) * 0.03) # estimate number of doublets. 3% multiplet rate for 3000-4000 cells recovered
infeT <- doubletFinder_v3(infeT, sct = TRUE, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30)
# [Infected] Doublet removal. Subset name comes from doubletFinder
infeT.sub <- subset(infeT, subset = DF.classifications_0.25_0.09_133 =="Singlet")

#***********************************
# Control-only analysis
#***********************************
#1. Dim reduction, clustering and differential expression
# SCTransform in place of normalization, finding variable features and scaling data.
ctrlT.sub <- SCTransform(ctrlT.sub, vars.to.regress = c("nFeature_RNA","percent.mt"))
ctrlT.sub <- RunPCA(ctrlT.sub, npcs = 50)
ctrlT.sub <- RunUMAP(ctrlT.sub, dims = 1:30)
# Clustering within control
ctrlT.sub <- FindNeighbors(ctrlT.sub, reduction = "pca", dims = 1:30)
ctrlT.sub <- FindClusters(ctrlT.sub, resolution = 0.6)
# Find markers for all
Idents(ctrlT.sub) <- "SCT_snn_res.0.6"
all.ctrl.markers <- FindAllMarkers(ctrlT.sub)

#2. T cell focused analysis
Idents(ctrlT.sub) <- "SCT_snn_res.0.6"
ctrlTonly <- subset(ctrlT.sub, idents = c("8","23","19","0","10","1","3"))
# Linear dim reduction
DefaultAssay(ctrlTonly) <- "SCT"
ctrlTonly <- RunPCA(ctrlTonly, npcs = 50)
# Non-linear dim reduction
ctrlTonly <- RunUMAP(ctrlTonly, reduction = "pca", dims = 1:20)
ctrlTonly <- FindNeighbors(ctrlTonly, reduction = "pca", dims = 1:20)
ctrlTonly <- FindClusters(ctrlTonly, resolution = 0.4)
# Find marker for all subtypes using res 0.4
DefaultAssay(ctrlTonly) <- "RNA"
Idents(ctrlTonly) <- "SCT_snn_res.0.4"
ctrlT.subtype.all.pos.markers <- FindAllMarkers(ctrlTonly, assay = "RNA", only.pos = TRUE)

# Combined Tc features for control only Tc subset
Tc.features <- c("nitr2b","nitr6b", # 0-Naive-1
                 "ccr7","cxcr4a", # 1-Naive-2
                 "mki67","top2a","stmn1a", # 6-cycling
                 "prf1.1","cxcr3.1","gzma", # 4-CTL-1
                 "nkl.2","ccr2","prf1.9", # 8-CTL-2
                 "il2rb","tnfsf14","ccl38.6", # 3-ccl38.6-high T cell
                 "foxp3a","ccr6a","ccl20a.3", # 5-Regulatory-like
                 "il11b","il4","il13", # 10-Th2/ILC2
                 "bcl11ba","rag1","rag2", # 9-rag1/2+ T
                 "ccl33.3","lcp2b", # 7-Lymphocyte-like
                 "cdh1","cldnb","cldnh"  # 2-Epithelial
                 )

# Rename clusters as cell types. Save ID as "Tc_identity"
DefaultAssay(ctrlTonly) <- "RNA"
DotPlot(ctrlTonly, features = Tc.features, group.by = "Tc_identity") + RotatedAxis()

# Remove epithelial cell cluster
Idents(ctrlTonly) <- "Tc_identity"
ctrlTonly.epirm <- subset(ctrlTonly, idents = "Epithelial cell", invert = TRUE)
Tc.features <- c("nitr2b","nitr6b", # 0-Naive-1
                 "ccr7","cxcr4a", # 1-Naive-2
                 "mki67","top2a","stmn1a", # 6-cycling
                 "prf1.1","cxcr3.1","gzma", # 4-CTL-1
                 "nkl.2","ccr2","prf1.9", # 8-CTL-2
                 "il2rb","tnfsf14","ccl38.6", # 3-ccl38.6-high T cell
                 "foxp3a","ccr6a","ccl20a.3", # 5-Regulatory-like
                 "il11b","il4","il13", # 10-Th2/ILC2
                 "bcl11ba","rag1","rag2", # 9-rag1/2+ T
                 "ccl33.3","lcp2b", # 7-Lymphocyte-like
                 )
DimPlot(ctrlTonly.epirm, group.by = "Tc_identity", order = rev(c("Naive T cell-like-1","Naive T cell-like-2","Cycling T cell","Cytotoxic T cell-1","Cytotoxic T cell-2","ccl38.6-high T cell", "Regulatory-like T cell","Th2 cell/ILC2" ,"rag1/2+ T cell", "Lymphocyte-like"))) # Figure S2G
DotPlot(ctrlTonly.epirm, features = Tc.features.epirm, group.by = "Tc_identity") + RotatedAxis() # Figure S2F
FeaturePlot(ctrlTonly.epirm, features = Tc.features.epirm) # Figure S2H

#3. Add detailed subtypes back to control
# Simplify Tc subset identities 
current.clid <- c("Epithelial cell","Naive T cell-like-2","Naive T cell-like-1","Cytotoxic T cell-2","Lymphocyte-like","Th2 cell/ILC2","Cycling T cell","rag1/2+ T cell","ccl38.6-high T cell","Cytotoxic T cell-1","Regulatory-like T cell")
new.clid <- c("Others-Epithelial cell", "T cell", "T cell", "T cell","Lymphocyte-like","T cell", "T cell", "T cell", "T cell", "T cell", "T cell")
ctrlTonly[["Tc_identity_simp"]] <- plyr::mapvalues(x = ctrlTonly$Tc_identity, from = current.clid, to = new.clid)
# Save cell name with each identity
Idents(ctrlTonly) <- "Tc_identity_simp"
Tc.Epi <- WhichCells(ctrlTonly, idents = "Others-Epithelial cell")
Tc.LL <- WhichCells(ctrlTonly, idents = "Lymphocyte-like")
Tc.Tc <- WhichCells(ctrlTonly, idents = "T cell")
# Renaming IDs from ctrlTonly to ctrlT.sub
ctrlT.sub <- SetIdent(ctrlT.sub, cells = Tc.Epi, value ="Others-Epithelial cell" )
ctrlT.sub <- SetIdent(ctrlT.sub, cells = Tc.LL, value ="Lymphocyte-like" )
ctrlT.sub <- SetIdent(ctrlT.sub, cells = Tc.Tc, value ="T cell" )
ctrlT.sub <- StashIdent(ctrlT.sub, save.name = "high_res_IDs")
# Visualization 
DimPlot(ctrlT.sub, group.by = "high_res_IDs", order = rev(c("T cell","Lymphocyte-like","Dendritic cell-like","B cell","Macrophage","Neutrophil","Erythrocyte","Thrombocyte","Superficial epithelial","Ionocyte","Intermediate epithelial","Basal epithelial","Mesenchymal","Lateral line-like","Epidermal mucous Cell","Others-Epithelial cell")),cols = DiscretePalette_scCustomize(num_colors = 16, palette = "ditto_seq")) # Figure 2J
# Major population marker check. simplified
all.major.ft.simp <- c(# "mki67","top2a", # cell cycle
                          "cd247l","lck", # general Tc
                          "ccl33.3", "lcp2b", # Lymphocyte like
                          "ctsbb","tlr7",# DC
                          "cd37","pax5", # B cell
                          "mpeg1.1","grn1", # Macrophage
                          "mpx","il6r", # Neutrophil
                          "hbba2","hemgn", # red blood cell
                          "thbs1b","fn1b", # Thrombocyte
                          "krt1-19d","cldne", # SE
                          "trpv6","foxi3b", # ionocyte 
                          "cldna","tp63", # IE
                          "cldn1","cldni", # BE
                          "vcana","clu",# Mes
                          "prox1a","prr15la", # lateral line
                          "agr2","cldnh" # Mucous
                          )
DefaultAssay(ctrlT.sub) <- "RNA"
ctrlT.sub[["relev_high_res_IDs"]] <- factor(ctrlT.sub$high_res_IDs, levels = rev(c("T cell","Lymphocyte-like","Dendritic cell-like","B cell","Macrophage","Neutrophil","Erythrocyte","Thrombocyte","Superficial epithelial","Ionocyte","Intermediate epithelial","Basal epithelial","Mesenchymal","Lateral line-like","Epidermal mucous Cell","Others-Epithelial cell")))
DotPlot(ctrlT.sub, features = all.major.ft.simp, group.by = "relev_high_res_IDs") + RotatedAxis() # Figure S2E
# Population structure
dittoBarPlot(object = ctrlT.sub, var = "high_res_IDs", group.by = "orig.ident", var.labels.reorder = c(15,9,3,1,10,12,5,16,14,7,6,2,11,8,4,13)) # Figure S2D
# All features on UMAP
FeaturePlot(ctrlT.sub, features = all.major.ft.simp) #

## Chemokine + receptor 
# Concise chemokine and receptors
receptors <- c("cxcr4a","cxcr4b","ccr7","cxcr5")
chemokines <- c("cxcl12a","cxcl12b","ccl19a.1","ccl19a.2","ccl19b","ccl25a","ccl25b","cxcl13")
DotPlot(ctrlT.sub, features = receptors ,group.by = "relev_high_res_IDs") + RotatedAxis() # Figure 2K
DotPlot(ctrlT.sub, features = chemokines ,group.by = "relev_high_res_IDs") + RotatedAxis() # Figure 2L


#**************************************
# Control-Infected integrative analysis
#**************************************
#1. Integrate with Seurat CCA
# normalization post doublet removal
ctrlT.sub <- SCTransform(ctrlT.sub, vars.to.regress = c("nFeature_RNA","percent.mt")) # if didn't run Control-only section
infeT.sub <- SCTransform(infeT.sub, vars.to.regress = c("nFeature_RNA","percent.mt"))

# Find features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(ctrlT.sub, infeT.sub))
# Find anchors between datasets
anchors <- FindIntegrationAnchors(object.list = list(ctrlT.sub, infeT.sub), anchor.features = features) #Using CCA (default) as the dim reduction approach
# Integrate
Tcells.combined <- IntegrateData(anchorset = anchors)

#2. Cell type identifications
# Set assay to integrated
DefaultAssay(Tcells.combined) <- "integrated"
# Perform linear dimensional reduction
Tcells.combined <- ScaleData(Tcells.combined)
Tcells.combined <- RunPCA(Tcells.combined, npcs = 50)
# Perform non-linear dimensional reduction
Tcells.combined <- RunUMAP(Tcells.combined, reduction = "pca", dims = 1:30)
Tcells.combined <- FindNeighbors(Tcells.combined, reduction = "pca", dims = 1:30)
Tcells.combined <- FindClusters(Tcells.combined, resolution = 0.5)
# Marker gene identifications
DefaultAssay(Tcells.combined) <- "RNA"
## Find all markers 
all.markers <- FindAllMarkers(Tcells.combined)

#3. T cell focused analysis
# Extract all T cells.
# | Cluster       | cd247l log2FC        | Pct.Exp  | lck log2FC | Pct.Exp |
# |---------------|---------------------:|---------:|-----------:|--------:|
# | 2             | 0.8                  | 0.33     | 0.6        | 0.4     |
# | 4             | -0.26                | 0.03     | 0.8        | 0.3     |
# | 6             | 0.7                  | 0.3      | 0.6        | 0.3     |
# | 9             | 0.4                  | 0.5      | 0.5        | 0.6     |
# | 15            | 1.2                  | 0.47     | NA         | NA      |
# | 16            | 1.0                  | 0.44     | NA         | NA      |
# | 17            | 1.1                  | 0.5      | 0.4        | 0.3     |
# Extract T cells from integrated data using log2FC > 0.3 and percent of expression >  0.3 for both cd247l and lck
Idents(Tcells.combined) <- "integrated_snn_res.0.5"
integrated.Tc <- subset(Tcells.combined, idents = c("2","6","9","17"))
# Split subset based on original identity (ctrl/infected)
tc.list <- SplitObject(integrated.Tc, split.by = "orig.ident")
# normalize and identify variable features for each dataset independently
tc.list <- lapply(X = tc.list, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("nFeature_RNA","percent.mt"))
})
# select features that are repeatedly variable across datasets for integration
inte.tc.ft <- SelectIntegrationFeatures(object.list = tc.list)
inte.tc.anchors <- FindIntegrationAnchors(object.list = tc.list, anchor.features = inte.tc.ft)
integrated.Tc.inte <- IntegrateData(anchorset = inte.tc.anchors)
# Perform integrated dim reduction
DefaultAssay(integrated.Tc.inte) <- "integrated"
integrated.Tc.inte <- ScaleData(integrated.Tc.inte)
integrated.Tc.inte <- RunPCA(integrated.Tc.inte, npcs = 50)
integrated.Tc.inte <- RunUMAP(integrated.Tc.inte, reduction = "pca", dims = 1:20)
integrated.Tc.inte <- FindNeighbors(integrated.Tc.inte, reduction = "pca", dims = 1:20)
integrated.Tc.inte <- FindClusters(integrated.Tc.inte, resolution = 0.4)
# Find markers in each group 
DefaultAssay(integrated.Tc.inte) <- "RNA"
Idents(integrated.Tc.inte) <- "integrated_snn_res.0.4"
all.Tc.markers <- FindAllMarkers(integrated.Tc.inte)
# Rename clusters as subset names and save in Identity
# Population structure within integrated stringent T cells
table(integrated.Tc.inte$Relev_ID, integrated.Tc.inte$orig.ident) # Figure 4J
# ｜Cell type          ｜Number of cells｜Ctrl｜Infected|
# ｜-------------------｜--------------:｜---:｜-------:|
# ｜Naïve T cell-like  ｜339            |217  |122     |
# ｜Cycling T cell     ｜299            |110  |189     |
# ｜Regulatory T cell-1｜335            |223  |112     |
# ｜Regulatory T cell-2｜182            |85   |97      |
# ｜Cytotoxic T cell-1 ｜155            |78   |77      |
# ｜Cytotoxic T cell-2 ｜129            |55   |74      |
# Markers for stringent T cell subsets only
Tc.stringent.ft <- c("ccr7","cxcr4a", # naive T
                     "mki67","top2a", # cycling T
                     "rorc", "il17ra1a", "foxp3a","ccr6a",# T reg 
                     "cd28l", # CTLA4
                     "cxcr3.1","gzma","prf1.1", # CTL-1
                     "nkl.2","ccr2","prf1.9", # CTL-2
                     "ifng1" # sign of activation
                     )
DimPlot(integrated.Tc.inte, reduction = "umap", group.by = "Relev_ID", cols = subtype.colors,order = rev(c("Naive T cell-like","Cycling T cell","Regulatory T cell-1","Regulatory T cell-2","Cytotoxic T cell-1","Cytotoxic T cell-2"))) # Figure S3D
DefaultAssay(integrated.Tc.inte) <- "RNA"
integrated.Tc.inte[["Relev_ID"]] <- factor(integrated.Tc.inte$Identity, levels = rev(c("Naive T cell-like","Cycling T cell","Regulatory T cell-1","Regulatory T cell-2","Cytotoxic T cell-1","Cytotoxic T cell-2")))
# Add orig.ident to the subtypes to show the control to infected comparison
integrated.Tc.inte[["Relev_ID_split"]] <- paste(integrated.Tc.inte$Relev_ID, integrated.Tc.inte$orig.ident, sep = "_")
# Relevel IDs
integrated.Tc.inte[["Relev_ID_split"]] <- factor(integrated.Tc.inte$Relev_ID_split, levels = rev(c("Naive T cell-like_ctrl","Naive T cell-like_infected","Cycling T cell_ctrl","Cycling T cell_infected","Regulatory T cell-1_ctrl","Regulatory T cell-1_infected","Regulatory T cell-2_ctrl","Regulatory T cell-2_infected","Cytotoxic T cell-1_ctrl","Cytotoxic T cell-1_infected","Cytotoxic T cell-2_ctrl","Cytotoxic T cell-2_infected")))
DotPlot(integrated.Tc.inte, features = Tc.stringent.ft, group.by = "Relev_ID_split") + RotatedAxis() # Figure S3E
FeaturePlot(integrated.Tc.inte, features = Tc.stringent.ft) #
# Naive T cell-like and cycling T cell comparisons 
Idents(integrated.Tc.inte) <- "Relev_ID"
VlnPlot(integrated.Tc.inte, features = "ccr7", idents = c("Naive T cell-like","Cycling T cell"), split.by = "orig.ident") # Figure 4K
VlnPlot(integrated.Tc.inte, features = "cd28l", idents = c("Naive T cell-like","Cycling T cell"), split.by = "orig.ident") # Figure 4L
VlnPlot(integrated.Tc.inte, features = "ifng1", idents = c("Naive T cell-like","Cycling T cell"), split.by = "orig.ident") # Figure 4M





