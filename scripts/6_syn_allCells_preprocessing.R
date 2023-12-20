#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
outName <- "allCells_syn"


################################
### BEGIN DATA PREPROCESSING ###
################################

### prepare healthy vs OA synovial fluid datasets

#load in 10x data and qc filter each sample
load10x(din = "./input_synovium/", dout = "./output/s1/", outName = "230929_rngr612_noMods", testQC = F, nFeature_RNA_high = 4000, nFeature_RNA_low = 100, percent.mt_high = 12.5, nCount_RNA_high = 25000, nCount_RNA_low = 500)

#integrate datasets into one object
seu.obj <- sctIntegrate(din = "./output/s1/synovium/", dout = "./output/s2/", outName = "230929_rngr612_noMods", vars.to.regress = "percent.mt", nfeatures = 2500)

#use the integrated object to identify ideal clustering resolution 
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230929_rngr612_noMods", test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete dimension reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230929_rngr612_noMods", final.dims = 50, final.res = 0.7, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.2, n.neighbors = 25, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "CAMP", "FLT3", "HLA-DRA", 
                                     "CD4", "MS4A1", "IL1B","CD68")
                       )

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


##############################
### END DATA PREPROCESSING ###
##############################

### NOTE: the output rds file: "./output/s3/230929_rngr612_noMods_res0.7_dims50_dist0.2_neigh25_S3.rds", contains the processed data that will be used of subset analysis of each major population. The below code is part of what was used to assign high level cell type identities to move forward with independent reclustering analysis

##################################
### BEGIN CELL TYPE ANNOTATION ###
##################################

#load in processed data & metadata
seu.obj <- readRDS("./output/s3/230929_rngr612_noMods_res0.7_dims50_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./syn_idents_09-29-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./syn_idents_09-29-2023.csv", groupBy = "clusterID", metaAdd = "intID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Normal","OA"))
outName <- "allCells_syn"

Idents(seu.obj) <- "majorID"
seu.obj <- RenameIdents(seu.obj, c("synoviocyte" = "Synoviocyte", "myeloid" = "Macrophage/DC", 
                                   "endothelial" = "Endothelial", "tcell" = "T cell",
                                  "myofibroblast" = "Myofibroblast", "fibroblast" = "Fibroblast")
                       )
seu.obj$majorID_pertyName <- Idents(seu.obj)


### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/viln/allCells/eqsf_n3_n3_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", "clusterID", "name", "cellSource"), skipEXPR = F,
               test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "FBLN1","ACAT2")
                          
                          )    


### Data supplemental - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "eqsyn_n1_n1",
                      outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)
                     


colArray <- read.csv("./syn_idents_09-29-2023.csv")
# colArray <- colArray %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colArray)*3)[ c( rep(FALSE, 2), TRUE ) ] )%>% arrange(clusterID)
# write.csv(colArray,"./syn_idents_09-29-2023.csv", row.names = F)


### Fig supp - plot inital cluster umap
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
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 1a: plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "intID",
#                 cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_majorID_UMAP.png", sep = ""), width = 7, height = 7)



### Fig XX: key feature plots
features <- c("PTPRC","CD3E","ESAM", 
              "DRA","AIF1","CD68",
              "COL14A1","IL13RA2", 
              "CD55","THY1", 
              "GPNMB","DEFB1",
              "SCGB1A1","COL28A1","ENSECAG00000036123","CXCL14")

colz <- c("#00BBDC","#00BBDC","#00BBDC",
          "#00BBDC","#00BBDC","#00BBDC",
          "#00C08B","#49B500","#49B500",
          "#49B500","#ED813E","#F8766D"
          )
title <- features
title <-  c("PTPRC (CD45)","CD3E","CTSW", 
              "DRA (MHCII)","AIF1 (Iba1)","CD68",
              "CSF3R","S100A12", 
              "FLT3","FCER1A", 
              "GPNMB","DEFB1",
              "COL1A2","MS4A1 (CD20)","FLT3")

p <- prettyFeats(seu.obj = seu.obj, pt.size = 0.00000001, nrow = 4, ncol = 4, title.size = 24, features = features, order = F, noLegend = T) 
ggsave(paste("./output/", outName, "/", outName, "_featPlots.png", sep = ""), width = 12, height = 12)


#######################################
### END INITAL CELL TYPE ANNOTATION ###
#######################################
