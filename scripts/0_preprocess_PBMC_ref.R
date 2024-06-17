#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

################################
### BEGIN DATA PREPROCESSING ###
################################

#load in 10x data and qc filter each sample
load10x(din = "./pbmc_data/", dout = "./output/s1_pbmc/", outName = "pbmc_data", testQC = F,
       nFeature_RNA_high = 3750, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 25000, nCount_RNA_low = 100, readCnts = T)

#integrate datasets into one object
seu.obj <- sctIntegrate(din = "./output/s1_pbmc/", dout = "./output/s2/", outName = "2023_09_11_pbmc_data", vars.to.regress = "percent.mt", nfeatures = 2500)

#use the integrated object to identify ideal clustering resolution 
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "2023_09_11_pbmc_data", test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete dimension reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "2023_09_11_pbmc_data", final.dims = 50, final.res = 0.4, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("CD3G", "TRAT1", "PRF1", "CTSW", 
                                     "MPEG1", "DRA", "CST3", "CD14", 
                                     "MS4A1", "CD79A", "LTC4S","TOP2A")
                       )

##############################
### END DATA PREPROCESSING ###
##############################

### NOTE: the output rds file: "./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds", contains the processed data that will be used of subset analysis of each major population. The below code is part of what was used to assign high level cell type identities to move forward with independent reclustering analysis

###########################################
### BEGIN CELL TYPE ANNOTATION TRANSFER ###
###########################################

#load in preprocesed data
seu.obj <- readRDS("./output/s3/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds")
outName <- "pbmc"

#split out sample IDs
namez <- lapply(colnames(seu.obj), function(x){
    strsplit(x, split = "\\.")[[1]][1]
})
unique(unlist(namez))

#load in cell type annotations obtained from Patel et al. -- https://github.com/BradRosenbergLab/equinepbmc/issues/3
cellID.df <- read.csv("./eq_idents.csv")
ct.list <- as.factor(cellID.df$clusterID_patel)
names(ct.list) <- cellID.df$eqName.barcode
seu.obj <- AddMetaData(seu.obj, metadata = ct.list, col.name = "clusterID_patel")


### Fig extra - plot annotations to ensure transfer looks to have occured correctly
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_patel",
           #    cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)

#load in additional cell type annotations
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./patel_idents.csv", groupBy = "clusterID_patel", metaAdd = "celltype.l1")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./patel_idents.csv", groupBy = "clusterID_patel", metaAdd = "celltype.l2")


### Fig extra - plot annotations to ensure transfer looks to have occured correctly
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "celltype.l1",
           #    cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = F
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_UMAP_ctl1.png", sep = ""), width = 7, height = 7)


### Fig extra - plot annotations to ensure transfer looks to have occured correctly
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "celltype.l2",
           #    cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = F
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_UMAP_ctl2.png", sep = ""), width = 7, height = 7)

#save the annotated Seurat object for future use
saveRDS(seu.obj, file = "./output/s3/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds")


#########################################
### END CELL TYPE ANNOTATION TRANSFER ###
#########################################
