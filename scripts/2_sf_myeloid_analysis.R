#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

#########################
# myeloid cell analysis
#########################
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
outName <- "myeloid"
datE <- "Aug_4_2023"
nfeatures <- 2500

#subset on t cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "myeloid")



table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- readRDS("./output/s2/Aug_4_2023_myeloid_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#during analysis several low qaulity clusters were identified and need to removed -- also remove a cluster of fibroblasts
seu.obj <- readRDS("./output/s3/Aug_8_2023_myeloid_2500_res1_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_myeloid_08-09-2023.csv", groupBy = "clusterID_sub", metaAdd = "keep")
outName <- "myeloid"
datE <- "Sep_1_2023"
nfeatures <- 2500


#subset on non-low quality cells
seu.obj <- subset(seu.obj,
                  subset = 
                  keep ==  "ok")

table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)
table(seu.obj$clusterID_sub)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = "40", algorithm = 3, prefix = "integrated_snn_res.")


seu.obj$clusterID_sub_orig <- seu.obj$clusterID_sub
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#subset on non-low quality cells
seu.obj$clusterID_sub_orig2 <- seu.obj$clusterID_sub
seu.obj <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "12")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)



###load in data
# seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_10_2023_myeloid_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
# seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_18_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Sep_1_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

Idents(seu.obj) <- "orig.ident"
seu.obj$cellSource <- ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal")

outName <- "myeloid"
datE <- "Sep_1_2023"

#set major Idents
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2" , "1" = "GPNMB_Macrophage", 
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



### Extra: Supplemental data: Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = paste(datE, outName, sep = "_"),
                      outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)


ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, dir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",datE,outName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID", "clusterID_sub", "majorID_sub", "name", "cellSource"), 
               skipEXPR = F, test = F,
                           feats = c("CD68", "FLT3", "CD1C")
                          
                          )

colArray <- read.csv("./sf_idents_08-01-2023.csv")
colArray.sub <- colArray[colArray$majorID == "myeloid",]
# colArray.sub$clusterID <- levels(seu.obj$clusterID_sub)
# colArray.sub$majorID <- levels(seu.obj$majorID_sub)
# write.csv(colArray.sub, "sf_idents_myeloid_12-05-2023.csv", row.names = F)

### Fig 2a: Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              cols = colArray.sub$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8)  + NoLegend() + theme(axis.title = element_blank(),
                                                                                  panel.border = element_blank())


axes <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
theme(axis.line = element_line(colour = "black", 
                               arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                             ends = "last", type = "closed"),
                              ),
      axis.title.y = element_text(colour = "black", size = 20),
      axis.title.x = element_text(colour = "black", size = 20),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
     )

p <- p + inset_element(axes,left= 0,
  bottom = 0,
  right = 0.25,
  top = 0.25,
                       align_to = "full")
ggsave(paste0("./output/", outName, "/", outName, "_rawUMAP_tiny.png"), width = 7, height = 7)


### Fig extra: Create UMAP by majorID
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "majorID_sub",
              pt.size = 0.25,
              cols = colArray.sub$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- formatUMAP(plot = pi) + NoLegend()
# p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = majorColors.df$labCol) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_majorSub_UMAP.png", sep = ""), width = 7, height = 7)


### Fig extra: Plot key feats
features <- c("PRF1",
              "GZMA", "GZMB",
              "SELL", "S100A12","IL1B",
              "GZMK",
              "CCL14", "C1QC",
              "MSR1","CSF1R","CCL3",
              "FLT3", "BATF3", "CADM1")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.25, title.size = 14, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 9, height = 15)


### Fig 2b: Create violin plots for key feats
features <- c("AIF1","DRA",
              "SELP", "CSF3R", "SELL","S100A12",
              "VCAN", "LYZ",
              "AQP9",
              "CCL14","MSR1","CSF1R",
              "CCL2","CCL8","LYVE1",
              "TPPP3",
              "MARCO","CD5L","GPNMB",
              "FLT3",
              "BATF3","DNASE1L3","IRF8",
              "CD1C","FCER1A",
              "TCF4","MS4A1",
              "LAMP3","CCR7",
              "TOP2A")

pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    group.by = "majorID_sub",
    combine = T,
    stack = T,
    cols = colArray.sub$newCol[c(8,9,7,5,4,3,2,6,1,10,12,11)],
    fill.by = "ident",
    flip = T,
    features = rev(features)
        ) & NoLegend() & theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_text(size = 7),
                               axis.text = element_text(size = 7),
                               axis.line = element_blank(),
                              strip.text.y.right = element_text(size = 7, hjust = 0)) + theme(plot.margin = margin(c(3,0,0,10)))

# plot <- prettyViln(plot = pi, colorData = NULL, nrow = length(features), ncol = 1)
ggsave(paste("./output/", outName, "/", outName, "selectViln.png", sep = ""), width = 3.75, height = 5)


### Fig 5e & Supp: M1 vs M2 signatures
#set modules
modulez <- list("Pro-inflammatory" = c("AZIN1", "CD38","CD86","CXCL10","FPR2","GPR18","IL12B","IL18","IRF5","NFKBIZ","NOS2","PTGS2","TLR4","TNF"),
                "Anti-inflammatory" = c("ALOX15", "ARG1", "CHIL3", "CHIL4","EGR2", "IL10","IRF4","KLF4","MRC1","MYC","SOCS2","TGM2")
                )

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = rev(features)
                    ) + coord_flip() + theme(plot.margin = margin(7, 7, 7, 7, "pt"),
                                                                   axis.title = element_blank(),
                                                                   legend.position = "right",
                                                                   legend.direction = "vertical",
                                                                   axis.text.x = element_text(angle=0, hjust = 0.5)
                                            ) + scale_y_discrete(position = "right") + scale_colour_continuous(name="Enrichment score", type = "viridis") + NoLegend()

ggsave(paste("./output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height = 1)


### Fig extra: plot each indiviudal feature
modulez <- c(list("Enrichment score" = names(modulez)), modulez)

plots <- lapply(modulez, function(x){
    majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = rev(unname(unlist(x)))
                    ) + coord_flip() + theme(axis.text.x = element_blank(),
                                                          axis.ticks.x = element_blank(),
                                                          legend.position = "right",
                                                                   legend.direction = "vertical",
                                                          axis.title = element_blank(),
                                                          plot.margin = margin(3, 0, 3, 0, "pt")
                                                         ) + scale_colour_viridis(option="magma", name='Average\nexpression', limits = c(-2.5,2.5))
    
})

    patch <- area()
    nrow <- length(modulez)
    ncol <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }


plots$`Enrichment score` <- plots$`Enrichment score` + theme(axis.text.x = element_text(angle=0, hjust = 0.5)
                            ) + scale_y_discrete(position = "right") + scale_colour_viridis() + guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))

p <- Reduce( `+`, plots ) +  plot_layout(guides = "collect", design = patch, 
                                                                             height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length)))))

ggsave(paste("./output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 6,height=7)


### Fig extra: umap by sample
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name",
              cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.5,
              label = FALSE,
              shuffle = TRUE
)
p <- formatUMAP(pi) + NoLegend() + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "umap_bySample.png", sep = ""), width =7, height = 7)


### Fig extra: stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "majorID_sub") +
scale_fill_manual(labels = levels(seu.obj$name), values = levels(seu.obj$colz)) + 
theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12)) 
# + scale_x_discrete(limits=c(),expand = c(0, 0))
ggsave(paste("./output/", outName,"/",outName, "_stackedBar.png", sep = ""), width =8, height = 5)


### Fig extra: Evlauate cell frequency by majorID_sub
freqy <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name", showPval = F,
               namez = "name", 
               colz = "colz"
              )

freqy <- freqy + ggpubr::stat_compare_means(method = "t.test",
                                   method.args = list(var.equal = F),
                                   aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)

ggsave(paste("./output/", outName, "/",outName, "_freqPlots_majorID_sub.png", sep = ""), width = 8.5, height = 9)


### Fig 2c: monte-carlo permitation test
library(scProportionTest)
log2FD_threshold <- 0.58

#method of downsampling profoundly impacts the intpretations -- rec down sample by name over cellSource
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "majorID_sub", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )
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
                                                         legend.position = "top",
                                                         plot.margin = margin(c(3,3,0,10))
                                                        ) + ylab("abundance change (log2FC)")


ggsave(paste("./output/", outName, "/",outName, "_propTest_majorID_sub-wDS.png", sep = ""), width = 3.5, height = 2, scale = 2 )


### Fig extra: pseudobulk analysis
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "OA", grepLabel = c("OA","Normal") #improve - fix this so it is more functional
                    )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^RPS", addLabs = NULL, mkDir = T
                     )


### Fig supp: DGE using wilcoxon
seu.obj$allCells <- "DGE analysis of myeloid cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("OA", "Normal"))
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = outName, labCutoff = 20, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )


### Fig extra: DGE using wilcoxon within each macrophage cell type
linDEG(seu.obj = seu.obj, groupBy = "majorID_sub", comparision = "cellSource", outDir = paste0("./output/", outName,"/linDEG/"), 
       outName = outName, labCutoff = 20, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )

### Fig 2d: run/plot gsea results
dgea.df <- read.csv("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/myeloid/myeloid_DGE analysis of myeloid cellsgeneList.csv")
geneListUp <- dgea.df %>% arrange(p_val_adj) %>% filter(p_val_adj < 0.01, avg_log2FC > 0) %>% pull(X)
geneListDwn <- dgea.df %>% arrange(p_val_adj) %>% filter(p_val_adj < 0.01, avg_log2FC < 0) %>% pull(X)

p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5", upCol = "red", dwnCol = "blue", upOnly = T, termsTOplot=35, trimTerm = T)+ theme(axis.title = element_text(size = 16))

pi <- p + scale_x_continuous(limits = c(-50,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "-log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))

ggsave(paste("./output/", outName, "/linDEG/", outName, "_all_enriched_terms_2.png", sep = ""), width = 7.5, height =9)



### Fig 2e: split feature plots highlighting key degs
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
features <- c("CCL2","S100A12","APOE","FABP5","SPP1","CCL5","MARCO")
p <- FeaturePlot(seu.obj.sub,features = features, pt.size = 0.000000001, split.by = "cellSource", order = T, by.col = F
                ) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                           axis.title.y.right = element_text(size = 16),
                                                           axis.ticks = element_blank(),
                                                           axis.title = element_blank(),
                                                           axis.line = element_blank(),
                                                           plot.title = element_text(size=16),
                                                           title = element_blank(),
                                                           plot.margin = unit(c(0, 0, 0, 0), "cm")
                                                          ) & scale_color_gradient(breaks = pretty_breaks(n = 3), 
                                                                                   limits = c(NA, NA), low = "lightgrey", 
                                                                                   high = "darkblue")  & scale_colour_viridis(option="magma", name='Expression')
ggsave(paste("./output/", outName, "/",outName, "_splitFeats.png", sep = ""), width = 10, height = 4)



