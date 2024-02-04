#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")


################################
### BEGIN DATA PREPROCESSING ###
################################

### prepare healthy vs OA synovial fluid datasets

#load in 10x data and qc filter each sample
load10x(din = "./input/", dout = "./output/s1/", outName = "230731_rngr612_noMods", testQC = F, nFeature_RNA_high = 4000, nFeature_RNA_low = 100, percent.mt_high = 12.5, nCount_RNA_high = 30000, nCount_RNA_low = 200)

#integrate datasets into one object
seu.obj <- sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = "230731_rngr612_noMods", vars.to.regress = "percent.mt", nfeatures = 2500)

#use the integrated object to identify ideal clustering resolution 
# seu.obj <- readRDS("./output/s2/230731_rngr612_noMods_seu.integrated.obj_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230731_rngr612_noMods", test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete dimension reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230731_rngr612_noMods", final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.2, n.neighbors = 25, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "CAMP", "FLT3", "HLA-DRA", 
                                     "CD4", "MS4A1", "IL1B","CD68")
                       )

#check QC parameters
outName <- "allCells"
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)

#load in metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj$orig.ident), "OA", "Normal")


### clus 7 looks bad - remove and repeat integration

#subset on all cells but cluster 7
seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID ==  "7")
table(seu.obj.sub$clusterID)
table(seu.obj.sub$orig.ident)

#reintegrate data
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230731_rngr612_noMods_2500", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

#test cluster resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230731_rngr612_noMods", test_dims = 45, algorithm = 3, prefix = "integrated_snn_res.")

#complete dimension reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230731_rngr612_noMods", final.dims = 45, final.res = 0.5, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.2, n.neighbors = 25, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "CAMP", "FLT3", "HLA-DRA", 
                                     "CD4", "MS4A1", "IL1B","CD68")
                       )

##############################
### END DATA PREPROCESSING ###
##############################

### NOTE: the output rds file: "./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds", contains the processed data that will be used of subset analysis of each major population. The below code is part of what was used to assign high level cell type identities to move forward with independent reclustering analysis

##################################
### BEGIN CELL TYPE ANNOTATION ###
##################################

#load in processed data and associated metadata 
#this includes cell type annotations that were manually determined after evaluting the block of code below
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Normal","OA"))
outName <- "preProcess"

Idents(seu.obj) <- "majorID"
seu.obj <- RenameIdents(seu.obj, c("tcell" = "T cell", "myeloid" = "Macrophage/DC", 
                                   "cycling" = "Cycling cells", "bcell" = "B cell")
                       )
seu.obj$majorID_pertyName <- Idents(seu.obj)


### Data supplemental - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "eqsf_n3_n3",
                      outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)


### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/viln/allCells/eqsf_n3_n3_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", "clusterID", "name", "cellSource"), skipEXPR = F,
               test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "FBLN1","ACAT2")
                          
                          )    
                     
#load in colors
colArray <- read.csv("./sf_idents_08-01-2023.csv")

### Fig extra - plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
               cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_preProcess_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig extra - plot umap by major ID
colArray.sub <- colArray[colArray$majCol == "yes",]
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
               group.by = "majorID_pertyName",
              cols = colArray.sub$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_preProcess_majorUMAP.png", sep = ""), width = 7, height = 7)


### Fig extra - plot key cell type defining features
features <- c("PTPRC","CD3E","CTSW", 
              "DRA","AIF1","CD68",
              "CSF3R","S100A12", 
              "FLT3","FCER1A", 
              "GPNMB","DEFB1",
              "COL1A2","MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj,pt.size = 0.00000001, nrow = 4, ncol = 4, title.size = 24, features = features, order = F, noLegend = T) 
ggsave(paste("./output/", outName, "/", outName, "_preProcess_featPlots.png", sep = ""), width = 12, height = 12, scale = 2)


#######################################
### END INITAL CELL TYPE ANNOTATION ###
#######################################
