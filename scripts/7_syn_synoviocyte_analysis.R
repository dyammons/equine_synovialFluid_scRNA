#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


### Load in synovium data
seu.obj <- readRDS("./output/s3/Nov_11_2023_syn_2500_res0.3_dims40_dist0.3_neigh30_S3.rds")
syn.id <- seu.obj$clusterID_sub



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


seu.obj <- AddMetaData(seu.obj, syn.id, col.name = "synID")

seu.obj$finalClusters <- ifelse(is.na(seu.obj$synID),
                               ifelse(seu.obj$majorID == "synoviocyte", "remove", as.character(seu.obj$majorID_pertyName)),
                               paste0("synoviocyte_",as.character(seu.obj$synID))
                               )


seu.obj <- subset(seu.obj, invert = T,
                 subset = finalClusters == "remove")

seu.obj.syn <- seu.obj



### Load in synovial fluid data pre subset on macrophage
seu.obj <- readRDS("./output/s3/Oct_30_2023_allCells_annotated.rds")
seu.obj.sf <- seu.obj


seu.list <- c(SplitObject(seu.obj.sf, split.by = "orig.ident"), SplitObject(seu.obj.syn, split.by = "orig.ident"))
table(seu.obj.sf$orig.ident)
table(seu.obj.syn$orig.ident)

outName <- "sf-syn"
datE <- "Nov_22_2023"
nfeatures <- 2500

#complete independent reclustering
seu.obj <- indReClus(seu.list = seu.list, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )


clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )



seu.obj <- readRDS("./output/s3/Nov_22_2023_sf-syn_2500_res1_dims40_dist0.3_neigh30_S3.rds")
outName <- "sf-syn"
datE <- "Nov_22_2023"


clusTrans.df <- table(seu.obj$clusterID_sub, seu.obj$celltype.l2) %>% melt() %>% group_by(Var.2) %>% top_n(., 1, value)

clusNames <- as.character(clusTrans.df$Var.2)
names(clusNames) <- clusTrans.df$Var.1

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, clusNames)
seu.obj$clusterID_sub_Trans <- Idents(seu.obj)

#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub_Trans", numOfFeats = 24, outName = paste0(datE,outName),
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

vilnPlots(seu.obj = seu.obj, resumeFile = "./output/viln/sf-syn/Nov_22_2023sf-syn_gene_list.csv", resume = T,  groupBy = "clusterID_sub_Trans", numOfFeats = 24, outName = paste0(datE,outName),
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )




