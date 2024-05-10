#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(scProportionTest)

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
seu.obj$colz <- factor(seu.obj$colz)


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
all.df <- read.csv("sf_idents_08-01-2023.csv")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
outName <- "allCells"

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


# Generate table with cell numbers for reviwer
as.data.frame(table(seu.obj$celltype.l2, seu.obj$name)) %>% 
    pivot_wider(names_from = "Var2", values_from = "Freq") %>%
    rename(`Cell type` = Var1) %>%
    write.csv(., "./output/supplementalData/cell_counts.csv", row.names = F)

# loop through dims for reviewer
#complete unsupervised clustering and dim reduction
lapply(c(10, 15, 20, 25, 30, 35, 40, 45, 50), function (x){
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "dim_test_for_rev", final.dims = x, 
                           final.res = 0.7, stashID = paste0("Dimensions_", x), return_obj = T, algorithm = 3,  
                           prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", 
                           saveRDS = F,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )    
    
    pi <- DimPlot(seu.obj, 
                  reduction = "umap", 
                  group.by = "celltype.l1",
                  pt.size = 0.25,
                  label = TRUE,
                  label.box = TRUE,
                  repel = T
     ) + ggtitle(paste0("Dimensions_", x))
    ggsave(paste("./output/", outName, "/", x, "_dim.png", sep = ""), width = 7, height = 7)
})


### Data supplemental 1 - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l1", numOfFeats = 24, outName = "supplemental_data_1", returnViln = F, 
          outDir = "./output/supplementalData/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
          assay = "RNA", min.pct = 0.25, only.pos = T)


### Data supplemental 2 - generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l2", numOfFeats = 24, outName = "supplemental_data_2", returnViln = F, 
          outDir = "./output/supplementalData/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
          assay = "RNA", min.pct = 0.25, only.pos = T)


### Export data to make supplemental table 4
df <- read.csv("./output/supplementalData/supplemental_data_4_gene_list.csv")
df <- df[ ,c("cluster","gene")] %>% group_by(cluster) %>% mutate(rowNum = row_number()) %>% filter(rowNum <= 25) %>% ungroup()
df <- spread(df, key = rowNum, value = gene)
write.csv(df, "./output/supplementalData/supplemental_table_base.csv", row.names = F)


#get colors for each cell type
colz.base <- c(tc.df$newCol, mye.df$newCol, 
               all.df[all.df$majorID == "cycling" & all.df$majCol == "yes", ]$newCol,
               all.df[all.df$majorID == "bcell" & all.df$majCol == "yes", ]$newCol)
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, c("Cycling_cells", "B_cell"))


### Fig extra - plot clusterID_final umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l2",
              cols = colz.base,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- formatUMAP(plot = pi, smallAxes = T)
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
 ) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.6, smallAxes = T) 
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
) + 
    labs(colour = "Cell source:") + 
    theme(
        legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_text(size = 14)
    ) + 
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 4), byrow = F))
fig1d <- formatUMAP(pi, smallAxes = T) 
ggsave(paste("./output/", outName, "/", outName, "_supp1a.png", sep = ""), width = 7, height = 7)


### Fig 1c - stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "celltype.l1") + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
               values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_fig1c.png", sep = ""), width =7, height = 5)



### Fig supp 1b: Evaluate cell frequency by cluster using monte carlo permutation
log2FD_threshold <- 0.58
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "celltype.l1", sample_1 = "Normal", sample_2 = "OA", sample_identity = "cellSource" )

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
ht <- daOR(seu.obj = seu.obj, groupBy = "celltype.l1", splitBy = "name", 
           outName = outName, outDir = paste0("./output/", outName, "/"), 
           t_map = TRUE, cluster_order = p$data$clusters)

png(file = paste0("./output/", outName, "/", outName, "_OR_heat.png"), width = 2000, height = 1000, res = 400)
par(mfcol = c(1, 1))         
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "left")
dev.off()


# Export supplemental data
table(seu.obj$celltype.l1, seu.obj$name) %>%
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
    write.csv(., "./output/supplementalData/allCells_permRes.csv", row.names = F)


### Fig supp 1c - Evlauate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 1, groupBy = "celltype.l1", comp = "cellSource", legTitle = "Cell source", refVal = "name", showPval = F,
              namez = "name", 
              colz = "colz"
              ) + NoLegend()

ggsave(paste("./output/", outName, "/",outName, "_supp3d.png", sep = ""), width = 12, height = 3)


### Fig 1d/e - DEG heatmap and scatter plot of DEGs in each major subset
seu.obj$celltype.l1 <- factor(gsub("_", " ", seu.obj$celltype.l1))
#complete dge analysis within each major subset
linDEG(seu.obj = seu.obj, groupBy = "celltype.l1", comparision = "cellSource", outDir = paste0("./output/", outName,"/linDEG/"), outName = "allCells", cluster = NULL, labCutoff = 10, noTitle = F, labsHide = "^ENSECAG", contrast = c("OA", "Normal"),
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, returnPlots = F
                  )


#load dge results and create heatmap of degs
files <- list.files(path = "./output/allCells/linDEG/", pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv, header = T)
df <- bind_rows(df.list, .id = "column_label")
df2 <- do.call(rbind, df.list) %>% mutate(contrast = "OA_vs_Normal")
colnames(df2)[1] <- "gene"
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


### Complete pseudobulk DGE by celltype.l1
createPB(seu.obj = seu.obj, groupBy = "celltype.l1", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("./output/", outName, "/pseudoBulk/"), min.cell = 5,
         grepTerm = "Normal", grepLabel = c("Normal", "OA")
)

df <- read.csv(paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"), row.names = 1)
df$horse <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
df <- df %>% mutate(horse = ifelse(sampleID == "Normal_3", "horse_4", paste0("horse_",horse)))
write.csv(df, paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("./output/", outName, "/pseudoBulk/celltype.l1_deg_metaData.csv"),
          padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          paired = T, pairBy = "horse", strict_lfc = F,
          idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


### Fig 1f: heatmap of dge results by major cell types
files <- lapply(levels(seu.obj$celltype.l1), function(x){
    list.files(path = paste0("./output/allCells/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

files <- unlist(files)
df.list <- lapply(files, read.csv, header = T)

write.csv(do.call(rbind, df.list), "./output/supplementalData/supplemental_data_3_gene_list.csv", row.names = F)

cnts_mat <- do.call(rbind, df.list)  %>% 
    mutate(
        direction = case_when(
            log2FoldChange <= -1 ~ "Down",
            log2FoldChange >= 1 ~ "Up",
            log2FoldChange < 1 & log2FoldChange > -1 ~ "n.s."
        )
    ) %>% 
    filter(!direction == "n.s.") %>%
    group_by(gs_base, direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = gs_base, values_from = nRow) %>% 
    as.matrix() %>% t()

colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"

#order by number of total # of DEGs
cnts_mat[is.na(cnts_mat)] <- 0
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        
cnts_mat <- rbind(cnts_mat, c(0, 0), c(0, 0), c(0, 0))
rownames(cnts_mat) <- NULL
rownames(cnts_mat) <- c("Cycling cells", "Dendritic cells", "Neutrophils", "Macrophage", 
                        "CD4 T cells", "CD8 T cells", "gd T cells", "B cells")
cnts_mat <- cnts_mat[ , c(2,1)]

png(file = paste0("./output/", outName, "/",outName, "_fig1f.png"), width=1250, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col = viridis(100), #circlize::colorRamp2(c(0,max(cnts_mat)), colors = c("white","red")),
              cluster_columns = F,
              column_title = gt_render(
                  paste0("<span style='font-size:14pt; color:black'>**# of DEGs**</span><br>",
                         "<span style='font-size:12pt; color:black'>(OA vs Normal)</span>")
              ),
              show_column_names = TRUE,
              column_title_side = "top",
              column_names_gp = gpar(fontsize = 14, col = "black"),
              column_names_rot = 0,
              column_names_centered = TRUE,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                  if(cnts_mat[i, j] > 10) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
                  } else if(cnts_mat[i, j] <= 10) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "grey80"))
                  }                          
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"), show_heatmap_legend = FALSE)
dev.off()


### Fig 1f: pct of cycling cells
seu.obj.sub <- subset(seu.obj,
                  subset = 
                  celltype.l1 ==  "Cycling_cells")

#spin off for Reviewer
#complete independent reclustering
seu.obj.sub <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "cycling_cell_subset", preSub = T, nfeatures = 2000,
                         vars.to.regress = "percent.mt"
                       )

#dim reduction and unsupervised clustering
seu.obj.sub <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "./output/s3/", outName = "cycling_cell_subset", final.dims = 30, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              repel = F
)
p1 <- formatUMAP(plot = pi) 
    
features <- c("CD68", "CD3E", 
              "FLT3", "IL23R")
p2 <- prettyFeats(seu.obj = seu.obj.sub, nrow = 2, ncol = 2, features = features, 
             color = "black", order = F, pt.size = 0.0000001, title.size = 18)
p <- p1 + p2
ggsave(paste("./output/", outName, "/", outName, "_extrA.png", sep = ""), width = 12, height = 6)


Idents(seu.obj.sub) <- "clusterID_sub"
seu.obj.sub <- RenameIdents(seu.obj.sub, c("0" = "tcell" , "1" = "myeloid", 
                                   "2" = "tcell", "3" = "tcell",
                                   "4" = "myeloid", "5" = "myeloid",
                                   "6" = "myeloid", "7" = "myeloid"
                       ))
seu.obj.sub$majorID_sub <- Idents(seu.obj.sub)

table(seu.obj.sub$majorID_sub, seu.obj.sub$cellSource) %>%
    as.data.frame() %>% 
    group_by(Var2) %>%
    mutate(prop = Freq / sum(Freq),
           tot = sum(Freq))

seu.obj.sub$allCells <- "All cells"
seu.obj.sub$allCells <- as.factor(seu.obj.sub$allCells)
createPB(seu.obj = seu.obj.sub, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("./output/", outName, "/pseudoBulk/"), min.cell = 5,
         grepTerm = "Normal", grepLabel = c("Normal", "OA")
)

df <- read.csv(paste0("./output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"), row.names = 1)
df$horse <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
df <- df %>% mutate(horse = ifelse(sampleID == "Normal_3", "horse_4", paste0("horse_",horse)))
write.csv(df, paste0("./output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"))

p <- pseudoDEG(metaPWD = paste0("./output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
          padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName, "/pseudoBulk/"), 
          outName = outName, returnVolc = T,
          paired = T, pairBy = "horse", strict_lfc = F,
          idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

pi  <- prettyVolc(plot = p[[1]], rightLab = "Up in OA", leftLab = "Up in Normal", arrowz = T, lfcCut = 1) + labs(x = "log2(Fold change)") + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank(),
                                                                                                                                  plot.title = element_text(face = "bold", hjust = 0.5, size = 20, vjust = 2)) + ggtitle("Cycling cells (OA vs Normal)")
ggsave(paste("./output/", outName, "/", outName, "_fig1e.png", sep = ""), width = 7, height = 7)


modulez <- list(
    CYCLING_DOWN = pi$data %>% filter(threshold == "Down") %>% pull(gene),
    CYCLING_UP = pi$data %>% filter(threshold == "Up") %>% pull(gene)
)

seu.obj.sub <- subset(seu.obj,
                      invert = T, 
                  subset = 
                  celltype.l1 ==  "Cycling_cells")

#complete module scoring
seu.obj.sub <- AddModuleScore(seu.obj.sub,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj.sub@meta.data)[grep("_score", names(seu.obj.sub@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
Idents(seu.obj.sub) <- "celltype.l1"
seu.obj.sub <- RenameIdents(seu.obj.sub, c("Macrophage" = "Macrophage" , "DC" = "Dendritic cells", 
                                   "Neutrophil" = "Neutrophils", "CD4" = "CD4 T cells",
                                   "gdT" = "gd T cells", "CD8" = "CD8 T cells",
                                   "Cycling_cells" = "Cycling cells", "B_cell" = "B cells"
                       ))
seu.obj.sub$celltype.l1_perty <- Idents(seu.obj.sub)
seu.obj.sub$celltype.l1_perty <- factor(seu.obj.sub$celltype.l1_perty, levels = rev(c("Dendritic cells", "Neutrophils", "Macrophage", 
                                                                          "CD4 T cells", "CD8 T cells", "gd T cells", 
                                                                          "B cells")))
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj.sub, groupBy = "celltype.l1_perty", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_celltypes.png", sep = ""), width = 4, height = 5)




#load dge results and create heatmap of degs
files <- c("./output/tcell/IL23R_gd_T2_vs_IL23R_gd_T1_all_genes.csv", 
           "./output/tcell/IL23R_gd_T1_vs_other_T_cells_all_genes.csv",
           "./output/tcell/IL23R_gd_T2_vs_other_T_cells_all_genes.csv",
           "./output/tcell/IL23R_gd_vs_other_T_cells_all_genes.csv",
           "./output/myeloid/cDC2_vs_cDC1_all_genes.csv"
          )

df.list <- lapply(files, read.csv, header = T)
write.csv(do.call(rbind, df.list), "./output/supplementalData/supplemental_data_5_gene_list.csv", row.names = F)


##################################
### END ALL CELL DATA ANALYSIS ###
##################################