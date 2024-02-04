#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

#################################################
### BEGIN LABEL TRANSFER TO ALL CELLS DATASET ###
#################################################


### Load in all cells and extract barcodes for Cycling_cells and B_cell
seu.obj.all <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")

Idents(seu.obj.all) <- "majorID"
seu.obj.all <- RenameIdents(seu.obj.all, c("tcell" = "T_cell", "myeloid" = "Macrophage-DC", 
                                   "cycling" = "Cycling_cells", "bcell" = "B_cell")
                       )
seu.obj.all$majorID_pertyName <- Idents(seu.obj.all)

seu.obj.all <- subset(seu.obj.all, subset = majorID_pertyName == "Cycling_cells" | majorID_pertyName == "B_cell")
seu.obj.all$majorID_pertyName <- droplevels(seu.obj.all$majorID_pertyName)


### Load in T cells and extract barcodes for each subset
seu.obj.tcell <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_10_2023_tcell_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
seu.obj.tcell <- loadMeta(seu.obj = seu.obj.tcell, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l2")
seu.obj.tcell <- loadMeta(seu.obj = seu.obj.tcell, metaFile = "./sf_idents_tcell_08-10-2023.csv", groupBy = "clusterID_sub", metaAdd = "celltype.l1")


### Load in myeloid cells and extract barcodes for each subset
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


### Load in all cells and transfer the cell type annotations

#load in all cells dataset
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
outName <- "allCells"

#remove unwanted cells
seu.obj <- subset(seu.obj, cells = c(colnames(seu.obj.myeloid), colnames(seu.obj.tcell), colnames(seu.obj.all)))

#transfer labels for celltype.l2 level
seu.obj <- AddMetaData(seu.obj, seu.obj.all$majorID_pertyName, col.name = "majorID")
seu.obj <- AddMetaData(seu.obj, metadata = as.factor(c(as.factor(seu.obj.myeloid$celltype.l2),as.factor(seu.obj.tcell$celltype.l2))), col.name = "celltype.l2")
seu.obj$celltype.l2 <- as.factor(ifelse(is.na(seu.obj$celltype.l2),as.character(seu.obj$majorID),as.character(seu.obj$celltype.l2)))
seu.obj$celltype.l2 %>% table()

#transfer labels for celltype.l1 level
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

#stash sample level metadata
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "colz")


### Complete independent reclustering on annotated & subset dataset
datE <- "Dec_19_2023"
nfeatures <- 2500
outName <- "allCells"

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check cluster resolution
# seu.obj <- readRDS("./output/s2/Aug_10_2023_tcell_2500_S2.rds")
# clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40"), algorithm = 3, prefix = "integrated_snn_res.")

#complete unsupervised clustering and dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.7, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#########################################
### END CELL TYPE ANNOTATION TRANSFER ###
#########################################

### NOTE: the output rds file: "./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds", contains the processed data that will be used for analysis of all cells.

####################################
### BEGIN ALL CELL DATA ANALYSIS ###
####################################

#load in the processed object
seu.obj <- readRDS("./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
tc.df <- read.csv("sf_idents_tcell_12-05-2023.csv")
mye.df <- read.csv("sf_idents_myeloid_12-05-2023.csv")

#fix annotation error
Idents(seu.obj) <- "celltype.l1"
seu.obj <- RenameIdents(seu.obj, c("Monocyte" = "Macrophage"))
seu.obj$celltype.l1 <- Idents(seu.obj)

#create new cluster ID number based on cluster size; smallest number (0) cooresponds to largest cluster
clusterID_major <- table(seu.obj$celltype.l1) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_major=row_number()-1) %>% arrange(clusterID_major) 

newID <- clusterID_major$clusterID_major
names(newID) <- clusterID_major$Var1
Idents(seu.obj) <- "celltype.l1"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_major <- Idents(seu.obj)

#get colors for each cell type
colz.base <- c(tc.df$newCol, mye.df$newCol, gg_color_hue(length(levels(seu.obj$celltype.l2)))[c(20,3)])
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, levels(seu.obj$celltype.l2)[c(2,6)])


### Fig extra - plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 1a - plot clusterID_major umap

#make colors for major cell types
colz.base <- gg_color_hue(length(unique(seu.obj$celltype.l1)))
colz.base <- colz.base[c(3,6,4,5,1,7,2,8)]
#create the plot
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_major",
              cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = F
 )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, smallAxes = T)
ggsave(paste("./output/", outName, "/", outName, "_fig1a.png", sep = ""), width = 7, height = 7)


### Fig 1b - key feature plots
features <- c("CD3E","AIF1","CD1C", 
              "CTSW","CD68","CSF3R",
              "GZMA","GPNMB","S100A12",
              "TRDC","MARCO","MS4A1",
              "KLRB1","FLT3","TOP2A")

colz.base <- gg_color_hue(length(levels(seu.obj$clusterID_major)))

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
ggsave(paste("./output/", outName, "/", outName, "_fig1b.png", sep = ""), width = 9, height = 15, scale = 1)


### Fig supp 1a - umap colorized by sample
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
ggsave(paste("./output/", outName, "/", outName, "_supp1a.png", sep = ""), width =7, height = 7)


### Fig 1c - stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "celltype.l1") + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_fig1c.png", sep = ""), width =7, height = 5)


### Fig 1d/e - stacked bar graph by colorID

#complete dge analysis within each major subset
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "celltype.l1", comparision = "cellSource", outDir = paste0("./output/", outName,"/linDEG/"), outName = "allCells", cluster = NULL, labCutoff = 10, noTitle = F, labsHide = "^ENSECAG", contrast = c("OA", "Normal"),
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, returnPlots = F
                  )


### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj, groupBy = "celltype.l1", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("./output/", outName, "/pseudoBulk/"), min.cell = 5,
         grepTerm = "Normal", grepLabel = c("Normal", "OA")
)

df <- read.csv(paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"), row.names = 1)
df$horse <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
df <- df %>% mutate(horse = ifelse(sampleID == "Normal_3", "horse_4", paste0("horse_",horse)))
write.csv(df, paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          paired = T, pairBy = "horse",
          idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

#load dge results and create heatmap of degs
files <- list.files(path = "/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/allCells/linDEG/", pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv, header = T)
df <- bind_rows(df.list, .id = "column_label")
df.res <- df %>% group_by(cellType) %>% summarize(cnts = n()) %>% as.data.frame()
df.res <- df.res[match(levels(seu.obj$celltype.l1), df.res$cellType), ]
df.res$cellType <- gsub("_", " ", df.res$cellType)
df.res$cnts <- as.numeric(df.res$cnts)

rownames(df.res) <- df.res$cellType
df.res$cellType <- NULL

mat.res <- df.res %>% arrange(desc(cnts)) %>% as.matrix()

png(file = paste0("./output/", outName, "/", outName, "deg_heatmap.png"), width=1000, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(as.matrix(mat.res),#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col=viridis(100),
              cluster_columns = F,
              column_title = "# of DEGs",
              show_column_names = FALSE,
              column_title_side = "top",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              cell_fun = function(j, i, x, y, width, height, fill) {
                  if(mat.res[i, j] > 100) {
                      grid.text(sprintf("%.0f", as.matrix(mat.res)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
                  } else if(mat.res[i, j] < 100) {
                      grid.text(sprintf("%.0f", as.matrix(mat.res)[i, j]), x, y, gp = gpar(fontsize = 14, col = "grey80"))
                  }                          
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()


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


##################################
### END ALL CELL DATA ANALYSIS ###
##################################