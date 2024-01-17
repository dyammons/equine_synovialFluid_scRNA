#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


### Load in synovium data -- subset on macrophage
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


outName <- "mac"
datE <- "Nov_22_2023"
nfeatures <- 2500

#subset on myeloid
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "myeloid")



table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

seu.obj.syn.mac <- seu.obj


### Load in synovial fluid data pre subset on macrophage
seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Sep_1_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

Idents(seu.obj) <- "orig.ident"
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal")


#set major Idents
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2" , "1" = "APOE_Macrophage", 
                                   "2" = "CD5L_Macrophage", "3" = "TPPP3_Macrophage",
                                   "4" = "CCL2_Macrophage", "5" = "cDC1",
                                   "6" = "Mo-Mac", "7" = "Neutrophil",
                                   "8" = "Monocyte", "9" = "pDC",
                                   "10" = "cycling_DC", "11" = "migDC"
                       ))


seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = levels(seu.obj$majorID_sub)[c(8,9,7,5,4,3,2,6,1,10,12,11)])


Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "DC" , "1" = "Macrophage", 
                                   "2" = "Macrophage", "3" = "Macrophage",
                                   "4" = "Macrophage", "5" = "DC",
                                   "6" = "Monocyte", "7" = "Neutrophil",
                                   "8" = "Neutrophil", "9" = "DC",
                                   "10" = "DC", "11" = "DC"
                       ))


seu.obj$majorID_sub2 <- Idents(seu.obj)


seu.obj.sf.mac <- seu.obj

seu.list <- c(SplitObject(seu.obj.sf.mac, split.by = "orig.ident"), SplitObject(seu.obj.syn.mac, split.by = "orig.ident"))
table(seu.obj.sf.mac$orig.ident)
table(seu.obj.syn.mac$orig.ident)


#complete independent reclustering
seu.obj <- indReClus(seu.list = seu.list, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )


clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

seu.obj <- readRDS("./output/s3/Nov_22_2023_mac_2500_res0.6_dims40_dist0.3_neigh30_S3.rds")


clusTrans.df <- table(seu.obj$clusterID_sub, seu.obj$majorID_sub) %>% melt() %>% group_by(Var.2) %>% top_n(., 1, value)

clusNames <- as.character(clusTrans.df$Var.2)
names(clusNames) <- clusTrans.df$Var.1

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, clusNames)
seu.obj$clusterID_sub_Trans <- Idents(seu.obj)



seu.obj$tissueSource <- ifelse(grepl("syn", seu.obj$orig.ident), "Synovium", "Synovial fluid")
table(seu.obj$clusterID_sub_Trans, seu.obj$tissueSource)


#                    Synovial fluid Synovium
#   APOE_Macrophage            2511      190
#   CCL2_Macrophage            1361      222
#   CD5L_Macrophage            2069      659
#   Mo-Mac                      704       18
#   Monocyte                    163        1
#   Neutrophil                  656       37
#   TPPP3_Macrophage           1563      171
#   cDC1                       1113       11
#   cDC2                       2692      176
#   cycling_DC                  124        6
#   migDC                        38        0
#   pDC                         179       79

table(seu.obj$clusterID_sub_Trans, seu.obj$orig.ident)

#                    syn_n1 syn_oa1
#   APOE_Macrophage     123      67
#   CCL2_Macrophage     137      85
#   CD5L_Macrophage     509     150
#   Mo-Mac               10       8
#   Monocyte              0       1
#   Neutrophil           12      25
#   TPPP3_Macrophage     98      73
#   cDC1                  2       9
#   cDC2                119      57
#   cycling_DC            6       0
#   migDC                 0       0
#   pDC                  34      45



### Fig extra - plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              #cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = T) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)



### Fig extra - plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "tissueSource",
              split.by = "tissueSource",
#               cols = c("mediumseagreen","mediumpurple1"),
              pt.size = 0.25,
              label = F,
              label.box = F,
              repel = TRUE
 )
p <- formatUMAP(plot = pi) + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            ) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_splitUMAP.png", sep = ""), width = 14, height = 7)



### Fig extra - plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub_Trans",
              split.by = "tissueSource",
#               cols = c("mediumseagreen","mediumpurple1"),
              pt.size = 0.25,
              label = F,
              label.box = F,
              repel = TRUE
 )
p <- formatUMAP(plot = pi) + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            )
ggsave(paste("./output/", outName, "/", outName, "_splitUMAP.png", sep = ""), width = 15, height = 7)

#explore the tissue specifity obeserved in cluster 8
seu.obj$split8 <- paste0(seu.obj$clusterID_sub, "_", seu.obj$tissueSource)
Idents(seu.obj) <- "split8"
FindMarkers(seu.obj, ident.1 = "8_Synovium", min.pct = 0.25, only.pos = T)


### Fig supp - stacked bar graph by majorID_sub
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "tissueSource", groupBy = "tissueSource", clusters = "clusterID_sub_Trans") +
#scale_fill_manual(labels = levels(seu.obj$name), values = levels(seu.obj$colz)) + 
theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12)) 
ggsave(paste("./output/", outName,"/",outName, "_stackedBar.png", sep = ""), width = 8, height = 5)



res.df <- as.data.frame(table(seu.obj$clusterID_sub_Trans, seu.obj$tissueSource))
res.df <- res.df[!res.df$Var2 == "Synovial fluid", ] %>% arrange(Freq)
res.df$Var1 <- factor(res.df$Var1, levels = res.df$Var1)

p <- ggplot(res.df, aes(x = Freq, y = Var1, fill = Var2)) +
            geom_bar(stat = "identity", width = 1, colour="white") +
            theme_classic() +
            theme(title = element_text(size= 14),
                  legend.title = element_blank(),
                  legend.text = element_text(size= 12),
                  legend.position = "top",
                  legend.direction = "horizontal",
                  plot.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key.size = unit(1,"line"),
                  axis.title.x = element_text(size = 18),
                  axis.title.y = element_blank(),
                  axis.text = element_text(size = 16),
                  plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
            ) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0,0)) + 
            xlab(label = "Number of cells") + 
            guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values="grey") + NoLegend()

ggsave(p, file = paste0("./output/", outName, "/", subName, "_suppfig5.png"), height = 7, width = 7)



outName <- "mac_syn"

seu.obj.sub <- subset(x = seu.obj, subset = tissueSource == "Synovium")

#method of downsampling profoundly impacts the intpretations -- rec dwn sample by name over cellSource
library(scProportionTest)

Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj.sub, downsample = min(table(seu.obj.sub$orig.ident)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "clusterID_sub_Trans", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )

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
