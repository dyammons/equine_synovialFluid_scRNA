#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(scProportionTest)

#################################################
### BEGIN LABEL TRANSFER TO ALL CELLS DATASET ###
#################################################

#load synoviocyte subset data
seu.obj <- readRDS("./output/s3/Nov_11_2023_syn_2500_res0.3_dims40_dist0.3_neigh30_S3.rds")
syn.id <- seu.obj$clusterID_sub

#load all cell synovium data & metadata
seu.obj <- readRDS("./output/s3/230929_rngr612_noMods_res0.7_dims50_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./syn_idents_09-29-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./syn_idents_09-29-2023.csv", groupBy = "clusterID", metaAdd = "intID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Normal","OA"))

Idents(seu.obj) <- "majorID"
seu.obj <- RenameIdents(seu.obj, c("synoviocyte" = "Synoviocyte", "myeloid" = "Macrophage/DC", 
                                   "endothelial" = "Endothelial", "tcell" = "T cell",
                                  "myofibroblast" = "Myofibroblast", "fibroblast" = "Fibroblast")
                       )
seu.obj$majorID_pertyName <- Idents(seu.obj)

#merge the metadata slots
seu.obj <- AddMetaData(seu.obj, syn.id, col.name = "synID")
seu.obj$finalClusters <- ifelse(is.na(seu.obj$synID),
                               ifelse(seu.obj$majorID == "synoviocyte", "remove", as.character(seu.obj$majorID_pertyName)),
                               paste0("synoviocyte_",as.character(seu.obj$synID))
                               )
seu.obj <- subset(seu.obj, invert = T,
                 subset = finalClusters == "remove")


#repeat unsupervised clustering and dim reduction (skipping reintegration here)
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "231220_rngr612_noMods", final.dims = 50, final.res = 0.7, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.2, n.neighbors = 25, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "CAMP", "FLT3", "HLA-DRA", 
                                     "CD4", "MS4A1", "IL1B","CD68")
                       )

#########################################
### END CELL TYPE ANNOTATION TRANSFER ###
#########################################

### NOTE: the output rds file: "./output/s3/231220_rngr612_noMods_res0.7_dims50_dist0.2_neigh25_S3.rds", contains the processed data that will be used of subset analysis of each major population.

####################################
### BEGIN ALL CELL DATA ANALYSIS ###
####################################

#load in processed data
seu.obj <- readRDS("./output/s3/231220_rngr612_noMods_res0.7_dims50_dist0.2_neigh25_S3.rds")
outName <- "allCells_syn"

#create new cluster ID number based on cluster size; smallest number (0) cooresponds to largest cluster
clusterID_major <- table(seu.obj$finalClusters) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_major=row_number()-1) %>% arrange(clusterID_major) 

newID <- clusterID_major$clusterID_major
names(newID) <- clusterID_major$Var1
Idents(seu.obj) <- "finalClusters"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_major <- Idents(seu.obj)

seu.obj$finalClusters <- ifelse(grepl("synoviocyte", seu.obj$finalClusters), paste0("Synoviocyte_c",seu.obj$clusterID_major), as.character(seu.obj$finalClusters))

### Data supplemental - generate violin plots of defining features

vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = "eqsyn_n1_n1",
                      outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", resume = T, resumeFile = "/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/viln/allCells/eqsyn_n1_n1_gene_list.csv",
                      min.pct = 0.25, only.pos = T)

vilnPlots(seu.obj = seu.obj, groupBy = "finalClusters", numOfFeats = 24, outName = "eqsyn_n1_n1",
                      outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)


### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, dir = "./output/cb_input/", 
               markers = "./output/viln/allCells/eqsyn_n1_n1_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", "clusterID", "name", "cellSource"), skipEXPR = F,
               test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "FBLN1","ACAT2")
                          
                          )

#set colors
colz.base <- c("#00C1A7","#00A6FF","#00BADE", "#64B200", "#00B5ED", "#00C0BB", "#619CFF", "#AEA200", "#DB8E00", "#B385FF", "#F8766D") 


### Fig 5a: plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_major",
              cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, smallAxes = T) & NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_fig5a.png", sep = ""), width = 7, height = 7)


### Fig 5b - key feature plots
features <- c("PTPRC","CD3E","AIF1","CD68",
              "CD34","ESAM","NOTCH3","ACTA2",
              "FBLN1","COL1A1","IGFBP2","CXCL14",
              "IL13RA2","DEFB1","DNASE1L3","CD55")

colz <- c("black",colz.base[10],colz.base[4],colz.base[4],
          colz.base[8],colz.base[8],colz.base[9],colz.base[9],
          colz.base[11],colz.base[11],colz.base[1],colz.base[1],
          colz.base[1],colz.base[1],colz.base[1],colz.base[1]
          )

title <-  c("PTPRC (CD34)","CD3E","AIF1 (Iba1)","CD68",
              "CD34","ESAM","NOTCH3","ACTA2",
              "FBLN1","COL1A1","IGFBP2","CXCL14",
              "IL13RA2","DEFB1","DNASE1L3","CD55")

p <- prettyFeats(seu.obj = seu.obj,pt.size = 0.00000001, nrow = 4, ncol = 4, title.size = 20, features = features, 
                 order = F, noLegend = T, titles = title, color = colz) 
ggsave(paste("./output/", outName, "/", outName, "_fig5b.png", sep = ""), width = 12, height = 12, scale = 1)


### Fig 5c: plot inital cluster umap
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))

log2FD_threshold <- 0.58

prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "finalClusters", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )
p <- permutation_plot(prop_test)  + theme(axis.title.y = element_blank(),
                                          legend.position = "top") + guides(colour = guide_legend("", nrow = 2, 
                                                                                                  byrow = TRUE)) + coord_flip()

res.df <- prop_test@results$permutation
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
                                               ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),
                                                         axis.title.x = element_blank(),
                                                         axis.title.y = element_text(size = 12),
                                                         legend.text = element_text(size = 12),
                                                         legend.title = element_text(size = 12),
                                                         legend.position = "top",
                                                         plot.margin = margin(c(3,3,0,30))
                                                        ) + ylab("abundance change (log2FC)")

ggsave(paste("./output/", outName, "/",outName, "_fig2c.png", sep = ""), width = 3.5, height = 2, scale = 2 )


