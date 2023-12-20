#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(scProportionTest)


############################################
### BEGIN SYNOVIOCYTE DATA PREPROCESSING ###
############################################

#load preprocessed allCells object & metadata
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

#set output params
outName <- "syn"
datE <- "Nov_11_2023"
nfeatures <- 2500

#subset on synoviocytes
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "synoviocyte")

#ensure subset was properly done
table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check ideal clustering resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

#dim reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.3, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

### Extra: Supplemental data: Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "eqsyn_n1_n1",
                      outDir = "./output/viln/syn/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)
                     

#plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
#                 cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


#suspected low quality cells identified, remove them and re-integrate
seu.obj <- subset(seu.obj, inver = T,
                  subset = 
                  clusterID_sub ==  "6")

#ensure subset was properly done
table(seu.obj$majorID)
table(seu.obj$clusterID_sub)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check cluster resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

#complete unsupervised clustering and dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.3, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


##########################################
### END SYNOVIOCYTE DATA PREPROCESSING ###
##########################################



#######################################
### BEGIN SYNOVIOCYTE DATA ANALYSIS ###
#######################################

seu.obj <- readRDS("./output/s3/Nov_11_2023_syn_2500_res0.3_dims40_dist0.3_neigh30_S3.rds")
syn.id <- seu.obj$clusterID_sub

### Extra: Supplemental data: Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "eqsyn_n1_n1",
                      outDir = "./output/viln/syn/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
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
                     

### 
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "clusterID_sub", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )

res.df <- prop_test@results$permutation
log2FD_threshold <- 0.58

res.df <- res.df %>% mutate(Significance = as.factor(ifelse(obs_log2FD < -log2FD_threshold & FDR < 0.01,"Down",
                                                            ifelse(obs_log2FD > log2FD_threshold & FDR < 0.01,"Up","n.s.")))
                           ) %>% arrange(obs_log2FD)

res.df$clusters <- factor(res.df$clusters, levels = c(res.df$clusters))


p <- ggplot(res.df, aes(x = clusters, y = obs_log2FD)) + 
geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
                    color = Significance)) + theme_bw() + geom_hline(yintercept = log2FD_threshold, 
                                                                     lty = 2) + geom_hline(yintercept = -log2FD_threshold, 
                                                                                           lty = 2) + 
geom_hline(yintercept = 0) + scale_color_manual(values = c("blue", "red","grey")
                                               ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                         axis.title.x = element_blank(),
                                                         legend.position = "top")


ggsave(paste("./output/", outName, "/",outName, "_propTest_clusterID-wDS.png", sep = ""), width = 6, height = 3)


#######################################
### END SYNOVIOCYTE DATA ANALYSIS ###
#######################################
