#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


#complete final all cells analysis
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
outName <- "allCells"


#########################

seu.obj.all <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")

Idents(seu.obj.all) <- "majorID"
seu.obj.all <- RenameIdents(seu.obj.all, c("tcell" = "T_cell", "myeloid" = "Macrophage-DC", 
                                   "cycling" = "Cycling_cells", "bcell" = "B_cell")
                       )

seu.obj.all$majorID_pertyName <- Idents(seu.obj.all)

seu.obj.all <- subset(seu.obj.all, subset = majorID_pertyName == "Cycling_cells" | majorID_pertyName == "B_cell")
seu.obj.all$majorID_pertyName <- droplevels(seu.obj.all$majorID_pertyName)


seu.obj.tcell <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_10_2023_tcell_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
seu.obj.tcell <- loadMeta(seu.obj = seu.obj.tcell, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l2")
seu.obj.tcell <- loadMeta(seu.obj = seu.obj.tcell, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l1")


seu.obj.myeloid <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Sep_1_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds")

Idents(seu.obj.myeloid) <- "clusterID_sub"
seu.obj.myeloid <- RenameIdents(seu.obj.myeloid, c("0" = "cDC2" , "1" = "GPNMB_Macrophage", 
                                   "2" = "CD5L_Macrophage", "3" = "TPPP3_Macrophage",
                                   "4" = "CCL2_Macrophage", "5" = "cDC1",
                                   "6" = "Mo-Mac", "7" = "Neutrophil",
                                   "8" = "Monocyte", "9" = "pDC",
                                   "10" = "cycling_DC", "11" = "migDC"
                       ))


seu.obj.myeloid$celltype.l2 <- Idents(seu.obj.myeloid)
seu.obj.myeloid$celltype.l2 <- factor(seu.obj.myeloid$celltype.l2, levels = levels(seu.obj.myeloid$celltype.l2)[c(8,9,7,5,4,3,2,6,1,10,12,11)])


Idents(seu.obj.myeloid) <- "clusterID_sub"
seu.obj.myeloid <- RenameIdents(seu.obj.myeloid, c("0" = "DC" , "1" = "Macrophage", 
                                   "2" = "Macrophage", "3" = "Macrophage",
                                   "4" = "Macrophage", "5" = "DC",
                                   "6" = "Monocyte", "7" = "Neutrophil",
                                   "8" = "Neutrophil", "9" = "DC",
                                   "10" = "DC", "11" = "DC"
                       ))


seu.obj.myeloid$celltype.l1 <- Idents(seu.obj.myeloid)




#remove unwanted cells
seu.obj <- subset(seu.obj, cells = c(colnames(seu.obj.myeloid), colnames(seu.obj.tcell), colnames(seu.obj.all)))


seu.obj <- AddMetaData(seu.obj, seu.obj.all$majorID_pertyName, col.name = "majorID")

seu.obj <- AddMetaData(seu.obj, metadata = as.factor(c(as.factor(seu.obj.myeloid$celltype.l2),as.factor(seu.obj.tcell$celltype.l2))), col.name = "celltype.l2")
seu.obj$celltype.l2 <- as.factor(ifelse(is.na(seu.obj$celltype.l2),as.character(seu.obj$majorID),as.character(seu.obj$celltype.l2)))

seu.obj$celltype.l2 %>% table()


seu.obj <- AddMetaData(seu.obj, metadata = as.factor(c(as.factor(seu.obj.myeloid$celltype.l1),as.factor(seu.obj.tcell$celltype.l1))), col.name = "celltype.l1")
seu.obj$celltype.l1 <- as.factor(ifelse(is.na(seu.obj$celltype.l1),as.character(seu.obj$majorID),as.character(seu.obj$celltype.l1)))

seu.obj$celltype.l1 %>% table()

#generate a cluster ID number in order of cluster size
clusterID_final <- table(seu.obj$celltype.l2) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

#stash the numerical ID
newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "celltype.l2"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)


seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")

# saveRDS(seu.obj,"./output/s3/Oct_30_2023_allCells_annotated.rds")
# saveRDS(seu.obj,"./output/s3/Nov_29_2023_allCells_annotated.rds")

#load in the annotated dataset
# seu.obj <- readRDS("./output/s3/Oct_30_2023_allCells_annotated.rds")
seu.obj <- readRDS("./output/s3/Nov_29_2023_allCells_annotated.rds")

tc.df <- read.csv("sf_idents_tcell_12-05-2023.csv")
mye.df <- read.csv("sf_idents_myeloid_12-05-2023.csv")

colz.base <- c(tc.df$newCol, mye.df$newCol, gg_color_hue(length(levels(seu.obj$celltype.l2)))[c(20,3)])
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, levels(seu.obj$celltype.l2)[c(2,6)])


### Fig supp: plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
#               cols = namez,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)

colz.base <- gg_color_hue(length(levels(seu.obj$celltype.l1)))
# colz.base <- gg_color_hue(length(levels(seu.obj$celltype.l2)))
### Fig supp: plot labeled umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l1",
              cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = F
 )
p <- formatUMAP(plot = pi)  + NoLegend() + theme(axis.title = element_blank(),
                                                 panel.border = element_blank(),
                                                 plot.margin = unit(c(-7, -7, -7, -7), "pt")
                                                )

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


ggsave(paste("./output/", outName, "/", outName, "_UMAP_ct.l1.png", sep = ""), width = 7, height = 7)


### Fig extra: plot inital cluster umap
pi <- DimPlot(seu.obj, 
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

# colArray.sub <- colArray[colArray$majCol == "yes",]
### Fig Fig 1a: umap by major ID
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
               group.by = "majorID_pertyName",
#               cols = colArray.sub$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_majorUMAP.png", sep = ""), width = 7, height = 7)


### Fig 1b: key feature plots
features <- c("CD3E","AIF1","CD1C", 
              "CTSW","CD68","CSF3R",
              "GZMA","GPNMB","S100A12",
              "TRDC","MARCO","MS4A1",
              "KLRB1","FLT3","TOP2A")

colz.base <- gg_color_hue(length(levels(seu.obj$majorID2)))

colz <- c("black","black",colz.base[4],
          colz.base[6],colz.base[3],colz.base[2],
          colz.base[6],colz.base[3],colz.base[2],
          colz.base[7],colz.base[3],colz.base[8],
          colz.base[7],colz.base[4],colz.base[1]
          )

title <-  c("CD3E","AIF1 (Iba1)","CD1C", 
              "CTSW","CD68","CSF3R",
              "GZMA","GPNMB","S100A12",
              "TRDC","MARCO","MS4A1 (CD20)",
              "KLRB1","FLT3","TOP2A")

p <- prettyFeats(seu.obj = seu.obj,pt.size = 0.00000001, nrow = 5, ncol = 3, title.size = 20, features = features, 
                 order = F, noLegend = T, titles = title, color = colz) 
ggsave(paste("./output/", outName, "/", outName, "_featPlots.png", sep = ""), width = 9, height = 15, scale = 1)




set.seed(12)
Idents(seu.obj) <- "name"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))

#create the plot
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "name",
              split.by = "name",
              pt.size = 0.25,
              cols = levels(seu.obj$colz),
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
)

pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 10.5, height = 7)


### Fig xx: key dot plot features
p <- majorDot(seu.obj = seu.obj, groupBy = "majorID2",
                  features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
                               "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
                               "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
                               "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_majorDot.png", sep = ""), width =8, height = 6)


### Fig 1d: umap by sample
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name",
              cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
fig1d <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_umap_bySample.png", sep = ""), width =7, height = 7)


### Fig 1e: stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "majorID2") + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_stackedBar.png", sep = ""), width =7, height = 5)


### Fig supp: stacked bar graph to supplement umap by sample
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID") +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1, byrow =T))
ggsave(paste("./output/", outName, "/", outName, "_supp_stackedBar.png", sep = ""), width =8, height = 12)


freqy <- freqPlots(seu.obj, method = 1, nrow= 2, comp = "cellSource", groupBy = "majorID2", legTitle = "Cell source",refVal = "name",
              namez = "name", 
              colz = "colz"
              ) + NoLegend()
ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 6, height = 4)


### Fig 1f: pct of cycling cells
seu.obj.sub <- subset(seu.obj,
                  subset = 
                  majorID ==  "cycling")


exp.df <- FetchData(object = seu.obj.sub, vars = c('cellSource', 'name', 'AIF1','CD3E'))

df <- exp.df %>% group_by(name) %>% mutate(aif_pos = ifelse(AIF1 > 0 & CD3E == 0, 1, 0),
                                     cd3e_pos = ifelse(CD3E > 0 & AIF1 == 0, 1, 0),
                                     total_cellz = n()
                                    ) %>% summarize(`AIF1+ cycling` = sum(aif_pos)/total_cellz,
                                                    `CD3E+ cycling` = sum(cd3e_pos)/total_cellz) %>% distinct() %>% as.data.frame() %>% melt(id.vars = "name") %>% mutate(cellSource = unlist(strsplit(as.character(name),"_"))[c(T,F)])

df$cellSource <- gsub("H", "Normal",df$cellSource)
df$name <- gsub("H", "Normal",df$name)

p <- ggplot(df, aes(y = value, x = cellSource)) + 
    labs(x = NULL, y = "% positive") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = NA, color = NA), 
          strip.text = element_text(face = "bold"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(color = "black")         )
    
    #note these are not corrected p-values
pi <- p + facet_wrap("variable", nrow = 1) + 
guides(fill = "none") +
geom_boxplot(aes_string(x = "cellSource"), alpha = 0.25, outlier.color = NA) + 
geom_point(size = 2, position = position_jitter(width = 0.25),
           aes_string(x = "cellSource", y = "value", color = "name")) +
labs(color = "Cell Source") +
ggpubr::stat_compare_means(method = "t.test",
                                   method.args = list(var.equal = F),
                                   aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3) + 
scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
theme(panel.grid.major = element_line(color = "grey", size = 0.25),
      text = element_text(size = 12)
          ) +                    
scale_color_manual(labels = levels(seu.obj$name), values = levels(seu.obj$colz))
                           
                                           
ggsave(paste("./output/", outName, "/", outName, "_pctOFcycling_facet.png", sep = ""), width =4.5, height = 3)


