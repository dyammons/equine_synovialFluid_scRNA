#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(scProportionTest)
library(UpSetR)

#######################################
### BEGIN T CELL DATA PREPROCESSING ###
#######################################

#load preprocessed allCells object
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
outName <- "tcell"
datE <- "Aug_8_2023"
nfeatures <- 2500

#subset on t cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "tcell")

#ensure subset was properly done
table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check ideal clustering resolution
# seu.obj <- readRDS("./output/s2/Aug_8_2023_tcell_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#dim reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

### A cluster of susspected cell doublets identified -- remove
seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_8_2023_tcell_2500_res1_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_tcell_08-08-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l2")
outName <- "tcell"
datE <- "Aug_10_2023"
nfeatures <- 2500

#remove susspected cluster of doublets
seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  celltype.l2 ==  "Doublets")
table(seu.obj.sub$celltype.l2)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check cluster resolution
# seu.obj <- readRDS("./output/s2/Aug_10_2023_tcell_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

#complete unsupervised clustering and dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.7, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#####################################
### END T CELL DATA PREPROCESSING ###
#####################################

#Note: the processed Seurat object is now stored as a .rds file, "./output/s3/Aug_10_2023_tcell_2500_res0.7_dims40_dist0.3_neigh30_S3.rds"

##################################
### BEGIN T CELL DATA ANALYSIS ###
##################################

#load in processed data and metadata
seu.obj <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_10_2023_tcell_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l1")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$celltype.l2 <- factor(seu.obj$celltype.l2, levels = levels(seu.obj$celltype.l2)[c(2,3,7,1,4,6,5,8,10,11,9,12)])
Idents(seu.obj) <- "orig.ident"
seu.obj$cellSource <- as.factor(ifelse(grepl("oa", seu.obj@meta.data$orig.ident), "OA", "Normal"))

#create new cluster ID number based on cluster size; smallest number (0) cooresponds to largest cluster
clusterID_final <- table(seu.obj$celltype.l2) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "celltype.l2"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)

#set file output specifications
outName <- "tcell"
datE <- "Aug_10_2023"


#summary data for cell type percentages in myeloid cells
table(seu.obj$celltype.l2, seu.obj$name) %>%
    as.data.frame() %>%
    separate(Var2, sep = "_", into = c("cellSource", NA),  remove = F) %>%
    group_by(Var2) %>%
    mutate(pct = prop.table(Freq)) %>%
    group_by(cellSource, Var1) %>%
    summarise(across(
        .cols = pct, 
        .fns = list(Mean = mean, MEDIAN = median, SD = sd, MIN = min, MAX = max), na.rm = TRUE, 
        .names = "{col}-{fn}"
    )) %>% 
    tidyr::pivot_longer(cols = where(is.double)) %>% 
    mutate(
        NAME = gsub("\\-.*", "", name),
        STAT = gsub(".*-", "", name),
        value = round((value * 100), 2)
    ) %>% 
    select(-name) %>% 
    pivot_wider(names_from = "NAME") %>% 
    as.data.frame() 


### Data supplemental - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l2", numOfFeats = 24, outName = "supplemental_data_7", returnViln = F, 
          outDir = "./output/supplementalData/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
          assay = "RNA", min.pct = 0.25, only.pos = T)

### Export data to make supplemental table 4
df <- read.csv("./output/supplementalData/supplemental_data_7_gene_list.csv")
df <- df[ ,c("cluster","gene")] %>% group_by(cluster) %>% mutate(rowNum = row_number()) %>% filter(rowNum <= 25) %>% ungroup()
df <- spread(df, key = rowNum, value = gene)
write.csv(df, "./output/supplementalData/supplemental_table_base.csv", row.names = F)

### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = "./output/supplementalData/supplemental_data_7_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                                   "clusterID", "clusterID_sub", "majorID", "celltype.l2", "name", "cellSource"), 
               skipEXPR = F, test = F,
               feats = c(
                   "CALY","CD40LG",
                   "SELL","LEF1","CCR7",
                   "FGL2","CTLA4","FOXP3",
                   "CTSW","GZMA", "CCL5","CCR5","IL2RB",
                   "GZMK","GZMM","IDO1",
                   "CRYBG2","CX3CR1","ZEB2",
                   "GNLY","KLRD1","KLRB1", 
                   "TRDC",
                   "SCART1","IL23R","KLRF1",
                   "BLK","CD160",
                   "IRF7","ISG20","ISG15"
               )
              )

#export the data for Zenodo
seu.obj <- cleanMeta(
    seu.obj = seu.obj, 
    metaSlot_keep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", 
        "nCount_SCT", "nFeature_SCT", "clusterID",
        "name", "cellSource", "celltype.l1", "celltype.l2"
    )
)
saveRDS(seu.obj, "./output/s3/eq_synovial_fluid_tcell_annotated.rds")

#load in colors
majorColors.df <- read.csv("sf_idents_tcell_08-10-2023.csv")
majorColors.df <- majorColors.df[!duplicated(majorColors.df$celltype.l2),]
majorColors.df <- majorColors.df %>% left_join(clusterID_final, by = c("celltype.l2" = "Var1"))
majorColors.df <- majorColors.df %>% arrange(celltype.l1) %>% mutate(newCol = gg_color_hue(nrow(majorColors.df)*8)[(nrow(majorColors.df)*8/2):(nrow(majorColors.df)*8)][ c( rep(FALSE, 3), TRUE )]) %>% arrange(clusterID_final)
# write.csv(majorColors.df, "sf_idents_tcell_12-05-2023.csv", row.names = F)
majorColors.df$labCol <- "black"


### Fig 3a: Create UMAP by clusterID_final
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
               cols = majorColors.df$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.6, smallAxes = T) & NoLegend()
ggsave(paste0("./output/", outName, "/", outName, "_fig3a.png"), width = 7, height = 7)


### Supp fig 3a: Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
#                cols = majorColors.df$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.6, smallAxes = T) & NoLegend()
ggsave(paste0("./output/", outName, "/", outName, "_supp3a.png"), width = 7, height = 7)


### Fig extra: Create UMAP by celltype.l2
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l2",
              pt.size = 0.25,
              cols = majorColors.df$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/",outName,"_ctl2_UMAP.png", sep = ""), width = 7, height = 7)


### Fig 3b: Plot key feats
features <-  c("CD3E","CCR7","LEF1","FOXP3",
              "GZMA", "TRDC", "SCART1", "IL23R")
p <- prettyFeats(seu.obj = seu.obj, pt.size = 0.00000001, nrow = 2, ncol = 4, title.size = 18, features = features, order = F, noLegend = T) 
ggsave(paste("./output/", outName, "/", outName, "_fig3b.png", sep = ""), width = 12, height = 6)


### Fig supp 3b: Plot key feats
features <-  c("CD4","CD8A")
p <- prettyFeats(seu.obj = seu.obj, pt.size = 0.00000001, nrow = 2, ncol = 1, title.size = 12, features = features, order = F, noLegend = T) 
ggsave(paste("./output/", outName, "/", outName, "_feat_CD4_CD8.png", sep = ""), width = 3, height = 6)


### Fig supp: plot enrichment scores
ecLists <- read.csv("../../k9_PBMC_scRNA/analysis/ecScores.csv", header = F)
modulez <- split(ecLists$V2, ecLists$V1)

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- paste0(names(modulez),"_SIG")

p <- majorDot(seu.obj = seu.obj, groupBy = "celltype.l2",
                     features = c("Naïve_SIG", "TH1_SIG","TH2_SIG", "TH17_SIG",
                                  "TREG_SIG")
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 24, "pt"),
                              legend.box="vertical") + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_supp_helper_ecScores.png", sep = ""), width = 6, height = 6)


### Fig extra: Create violin plots for key feats
features <- c("CTSW","TRAT1","TRDC", 
              "GZMA", "GZMM",
              "GZMK")
pi <- VlnPlot(object = seu.obj,
              pt.size = 0,
              same.y.lims = F,
              group.by = "clusterID_final",
              combine = T,
              cols = majorColors.df$newCol,
              stack = T,
              fill.by = "ident",
              flip = T,
              features = features
             ) + NoLegend() + theme(axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank()
                                   )
ggsave(paste("./output/", outName, "/", outName, "_selectViln.png", sep = ""), width = 5, height = 6)


### Fig extra: make dot plots of key features
features <- c("CALY","CD40LG",
              "SELL","LEF1","CCR7",
              "FGL2","CTLA4","FOXP3",
              "CTSW","GZMA", "CCL5","CCR5","IL2RB",
              "GZMK","GZMM","IDO1",
              "CRYBG2","CX3CR1","ZEB2",
              "GNLY","KLRD1","KLRB1", 
              "TRDC",
              "SCART1","IL23R","KLRF1",
              "BLK","CD160",
              "IRF7","ISG20","ISG15")
pi <- majorDot(seu.obj = seu.obj, groupBy = "celltype.l2",
               features = features
              ) + theme(axis.title = element_blank(),
                        axis.text = element_text(size = 12)
                       )
ggsave(paste("./output/", outName, "/", outName, "_majorDot.png", sep = ""), width = 10, height = 5)


### Fig supp 3b - reference map using eq PBMC data

#load in preprocessed and annotated equine pbmc data from patel et. al
reference <- readRDS("./output/s3/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds")
reference$celltype.l2 %>% table()
Idents(reference) <- "celltype.l2"
reference <- subset(reference,
                  idents = c("T CD4+ non-naïve","T CD4+ naïve","PRF1+ non-annotated","T CD8+ memory",
                            "T CD4+ cytotoxic","T CD8+ naïve","NK","T gd"))

reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", #reference.reduction = "umap",
    dims= 1:50 #dims= 1:2
)

predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l2,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              repel = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + theme(axis.title = element_blank(),
                             panel.border = element_blank()
                            ) + NoLegend() + coord_cartesian(clip = "off")
ggsave(paste("./output/", outName,"/",outName, "_supp3b.png", sep = ""), width = 7, height = 7)

#inspect numbers
table(seu.obj$predicted.id,seu.obj$celltype.l2)


### Fig supp 3c - dge analysis IL23R_gd_T2 vs IL23R_gd_T1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "IL23R_gd_T2", idents.2 = "IL23R_gd_T1", 
                      bioRep = "name",padj_cutoff = 0.01, lfcCut = 1, labSize = 4.5, strict_lfc = T, 
                      minCells = 5, outDir = paste0("./output/", outName, "/"), title = "IL23R_gd_T2_VS_IL23R_gd_T1", 
                      idents.1_NAME = "IL23R gd T2", idents.2_NAME = "IL23R gd T1", returnVolc = T, 
                      doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "IL23R gd T2", lfcCut = 1,
                 leftLab = "IL23R gd T1", arrowz = T) + labs(x = "log2(FC)") + theme(axis.title = element_text(size = 24),
                                                                                     axis.text = element_text(size = 18),
                                                                                     plot.title = element_blank(),
                                                                                     legend.position = c(0.10, 0.85)
                                                                                    ) + ggtitle("IL23R gd T2 vs IL23R gd T1")

ggsave(paste("./output/", outName, "/", outName, "_supp3c.png", sep = ""), width = 7, height = 7)


### Fig supp 3d - Evlauate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 3, groupBy = "celltype.l2", comp = "cellSource", legTitle = "Cell source", refVal = "name", showPval = F,
              namez = "name", 
              colz = "colz"
              ) + NoLegend()
ggsave(paste("./output/", outName, "/",outName, "_supp3d.png", sep = ""), width = 6, height = 5)


### Fig supp 3e - umap by sample
Idents(seu.obj) <- "orig.ident"
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
p <- formatUMAP(pi) + theme(legend.justification="center",
                            legend.position = "top"
                            ) + guides(color = guide_legend(nrow = 1, override.aes = list(size = 5)))
ggsave(paste("./output/", outName, "/", outName, "_supp3e.png", sep = ""), width = 7, height = 7)


### Fig 3c: Evaluate cell frequency by cluster using monte carlo permutation
log2FD_threshold <- 0.58
Idents(seu.obj) <- "name"
set.seed(12)

# seu.obj <- subset(x = seu.obj, subset = name %in% c("Normal_1", "OA_1"))
# seu.obj <- subset(x = seu.obj, subset = name %in% c("Normal_2", "OA_2"))
# seu.obj <- subset(x = seu.obj, subset = name %in% c("Normal_1", "Normal_2", "Normal_3", "OA_3"))
# Idents(seu.obj) <- "cellSource"
# seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

seu.obj$name <- droplevels(seu.obj$name)
Idents(seu.obj) <- "name"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "celltype.l2", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )

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
                                               ) + theme(axis.text.x = element_blank(),
                                                         axis.title.x = element_blank(),
                                                         legend.position = c(0.2, 0.9),
                                                         legend.title = element_blank(),
                                                         plot.margin = margin(c(7,7,21,7)),
                                                         legend.key = element_rect(fill = 'transparent', colour = NA),
                                                         legend.direction = 'vertical',
                                                         legend.background = element_rect(fill = 'transparent', colour = NA),
                                                         panel.background = element_rect(fill = 'transparent', colour = NA),
                                                         plot.background = element_rect(fill = "transparent", colour = NA)
                                                        ) + ylab(expression(atop(bold("abundance log2(FC)"),atop(italic("(OA versus Normal)"))))) + guides(color = guide_legend(nrow=1))

ggsave(paste("./output/", outName, "/",outName, "_fig3c.png", sep = ""), width = 3, height = 1, scale = 2 )


### Plot the assocaited heatmap
ht <- daOR(seu.obj = seu.obj, groupBy = "celltype.l2", splitBy = "name", 
           outName = outName, outDir = paste0("./output/", outName, "/"), 
           t_map = TRUE, cluster_order = p$data$clusters)

png(file = paste0("./output/", outName, "/", outName, "_OR_heat.png"), width = 2000, height = 1000, res = 400)
par(mfcol = c(1, 1))         
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "left")
dev.off()


#column_title_gp = gpar(fontsize = 14, fontface = "bold")
# Export supplemental data
table(seu.obj$celltype.l2, seu.obj$name) %>%
    as.data.frame() %>%
    separate(Var2, sep = "_", into = c("cellSource", NA),  remove = F) %>%
    group_by(Var2) %>%
    mutate(pct = round(prop.table(Freq), 4)) %>%
    group_by(cellSource, Var1) %>%
    summarise(across(
        .cols = pct, 
        .fns = list(MEAN = mean, SD = sd), na.rm = TRUE, 
        .names = "{col}-{fn}"
    )) %>% 
    tidyr::pivot_longer(cols = where(is.double)) %>% 
    mutate(
        NAME = gsub("\\-.*", "", name),
        STAT = gsub(".*-", "", name),
        value = round((value * 100), 2)
    ) %>% 
    
    select(-name, -NAME) %>% 
    pivot_wider(names_from = c(STAT, cellSource), values_from  = value) %>%
    left_join(p$data[ , -c(2, 3)], by = c("Var1" = "clusters")) %>%
    rename(`Cell Type` = Var1) %>%
    write.csv(., "./output/supplementalData/tcell_permRes.csv", row.names = F)


### Fig 3d - gene sig comparing IL23R_gd_T1 to other T cells
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "IL23R_gd_T1", idents.2 = NULL, bioRep = "name",padj_cutoff = 0.01, lfcCut = 1, topn = c(10, 10), labSize = 4,
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "IL23R_gd_T1_VS_otherT", idents.1_NAME = "IL23R gd T1", idents.2_NAME = "other T cells", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24, strict_lfc = T
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = NULL, lfcCut = 1,
                 leftLab = NULL, arrowz = F) + labs(x = "log2(FC)") + theme(axis.title = element_text(size = 24),
                                                                                                      axis.text = element_text(size = 18),
                                                                                                      plot.title = element_text(size = 28),
                                                                            legend.text = element_text(size= 16),
                                                                                     legend.position = c(0.15, 0.75)
                                                                                                      
                                                     ) + ggtitle("IL23R gd T1 vs other T cells") + NoLegend()

ggsave(paste("./output/", outName, "/", outName, "_fig3d.png", sep = ""), width = 7, height = 7)


### Fig extra - GSEA of gene sig comparing IL23R_gd_T1 to other T cells
degs.df <- read.csv("./output/tcell/IL23R_gd_T1_vs_other_T_cells_all_genes.csv")
degs.df %>% mutate(rection = ifelse(log2FoldChange > 0, "UP", "DWN")) %>% group_by(rection) %>% summarize(cnts = n())
up.genes <- degs.df %>% filter(log2FoldChange > 0) %>% pull(gene)
p <- plotGSEA(geneList = up.genes, species = "equine", geneListDwn = NULL, category = "C5", upCol = "red", dwnCol = "blue", upOnly = T, termsTOplot=1000, lolli = T, size = 3.5)
pi <- p + theme(axis.title=element_text(size = 16)) + scale_x_continuous(limits = c(-10,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Gene ontology (IL23R gd T1)") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste("./output/", outName,"/", outName, "_gd_il23r_T1_enriched_terms.png", sep = ""), width = 7, height =7)

gsea_res_1 <- pi$data %>% mutate(Contast = "IL23R_gd_T1_vs_other_T_cells")


### Fig 3e - gene sig comparing IL23R_gd_T2 to other T cells
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "IL23R_gd_T2", idents.2 = NULL, bioRep = "name",padj_cutoff = 0.01, lfcCut = 1, labSize = 4, topn = c(10,10),
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "IL23R_gd_T2_VS_otherT", idents.1_NAME = "IL23R gd T2", idents.2_NAME = "other T cells", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24, strict_lfc = T
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = NULL, lfcCut = 1,
                 leftLab = NULL, arrowz = F) + labs(x = "log2(FC)") + theme(axis.title = element_text(size = 24),
                                                                                                      axis.text = element_text(size = 18),
                                                                                                      plot.title = element_text(size = 28),
                                                                            legend.text = element_text(size= 16),
                                                                                     legend.position = c(0.15, 0.9)
                                                                                                      
                                                     ) + ggtitle("IL23R gd T2 vs other T cells") + NoLegend()

ggsave(paste("./output/", outName, "/", outName, "_fig3e.png", sep = ""), width = 7, height = 7)


### Fig extra - GSEA of gene sig comparing IL23R_gd_T2 to other T cells
degs.df <- read.csv("./output/tcell/IL23R_gd_T2_vs_other_T_cells_all_genes.csv")
degs.df %>% mutate(rection = ifelse(log2FoldChange > 0, "UP", "DWN")) %>% group_by(rection) %>% summarize(cnts = n())
up.genes <- degs.df %>% filter(log2FoldChange > 0) %>% pull(gene)
p <- plotGSEA(geneList = up.genes, geneListDwn = NULL, species = "equine", category = "C5", upCol = "red", dwnCol = "blue", upOnly = T, termsTOplot=1000, lolli = T, size = 3.5)
pi <- p + theme(axis.title=element_text(size = 16)) + scale_x_continuous(limits = c(-18,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Gene ontology (IL23R gd T2)") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste("./output/", outName,"/", outName, "_gd_il23r_T2_enriched_terms.png", sep = ""), width = 7, height = 7)

gsea_res_2 <- pi$data %>% mutate(Contast = "IL23R_gd_T2_vs_other_T_cells")


### Fig supp 3f - upset plot of gd T cell gene signature overlap
#load in dge results
IL23R_gd_T1.df <- read.csv("./output/tcell/IL23R_gd_T1_vs_other_T_cells_all_genes.csv")
IL23R_gd_T2.df <- read.csv("./output/tcell/IL23R_gd_T2_vs_other_T_cells_all_genes.csv")

#organize the data
upSet.df <- as.data.frame(unique(c(IL23R_gd_T1.df$gene,IL23R_gd_T2.df$gene)))
colnames(upSet.df) <- "gene"
upSet.df$`IL23R gd T1` <- as.integer(ifelse(upSet.df$gene %in% IL23R_gd_T1.df$gene, 1, 0))
upSet.df$`IL23R gd T2` <- as.integer(ifelse(upSet.df$gene %in% IL23R_gd_T2.df$gene, 1, 0))

#plot the upset
png(file = paste0("./output/", outName, "/", outName, "_supp3f.png"), width=1500, height=2000, res=400)
par(mfcol=c(1,1))     
p <- upset(upSet.df, sets = colnames(upSet.df)[2:3], mb.ratio = c(0.65, 0.35))
p
dev.off()


### Fig extra - dge analysis IL23R_gd_T1 & IL23R_gd_T2 vs other T cells
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = c("IL23R_gd_T1", "IL23R_gd_T2"), idents.2 = NULL, bioRep = "name",
                      padj_cutoff = 0.01, lfcCut = 1, topn = c(10,10), labSize = 4,
                      minCells = 5, outDir = paste0("./output/", outName, "/"), 
                      title = "IL23R_gd_VS_otherT", idents.1_NAME = "IL23R gd", idents.2_NAME = "other T cells", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24, strict_lfc = T
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = NULL, lfcCut = 1,
                 leftLab = NULL, arrowz = F) + labs(x = "log2(FC)") + theme(axis.title = element_text(size = 24),
                                                                                                      axis.text = element_text(size = 18),
                                                                                                      plot.title = element_text(size = 28),
                                                                            legend.text = element_text(size= 16),
                                                                                     legend.position = c(0.15, 0.75)
                                                                                                      
                                                     ) + ggtitle("IL23R gd vs other T cells")

ggsave(paste("./output/", outName, "/", outName, "_il23r_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig 3f - GSEA of dge results IL23R_gd_T1 & IL23R_gd_T2 vs other T cells
degs.df <- read.csv("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/tcell/IL23R_gd_vs_other_T_cells_all_genes.csv")
degs.df %>% mutate(rection = ifelse(log2FoldChange > 0, "UP", "DWN")) %>% group_by(rection) %>% summarize(cnts = n())
up.genes <- degs.df %>% filter(log2FoldChange > 0) %>% pull(gene)
p <- plotGSEA(geneList = up.genes, species = "equine", geneListDwn = NULL, category = "C5", upCol = "red", dwnCol = "blue", upOnly = T, lolli = T,  termsTOplot = 20, size = 3.75, saveRes = paste("./output/", outName,"/", outName, "_fig3f_gsea_source.csv", sep = ""))

pi <- p + theme(axis.title = element_text(size = 24),
                axis.text = element_text(size = 18)
               ) + scale_x_continuous(limits = c(-17,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)),name = "-log10(p.adj)") + labs(title = "Gene ontology", subtitle = "(IL23R gd T1/2 vs other T cells)") + theme(plot.title = element_text(size = 28, hjust = 0.5),  plot.subtitle = element_text(size = 18, hjust = 0.5))

ggsave(paste("./output/", outName,"/", outName, "_fig3f.png", sep = ""), width = 8, height = 7)

gsea_res_3 <- pi$data %>% mutate(Contast = "IL23R_gd_T1_AND_T2_vs_other_T_cells")

gsea_res <- rbind(gsea_res_1, gsea_res_2, gsea_res_3)
write.csv(gsea_res, "./output/supplementalData/supplemental_data_8_gene_list.csv", row.names = F)

################################
### END T CELL DATA ANALYSIS ###
################################