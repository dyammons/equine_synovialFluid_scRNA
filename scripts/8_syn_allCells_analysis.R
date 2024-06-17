#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
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
colz.base <- c("#00C1A7","#00A6FF","#00BADE", "#64B200", "#00B5ED", "#00C0BB", "#619CFF", "#AEA200", "#DB8E00", "#B385FF", "#F8766D") 
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


#summary data for cell type percentages in myeloid cells
table(seu.obj$finalClusters, seu.obj$name) %>%
    as.data.frame() %>%
    group_by(Var2) %>%
    mutate(pct = prop.table(Freq)) %>%
    group_by(Var1) %>%
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
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = "supplemental_data_8", returnViln = F, 
          outDir = "./output/supplementalData/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
          assay = "RNA", min.pct = 0.25, only.pos = T)

vilnPlots(seu.obj = seu.obj, groupBy = "finalClusters", numOfFeats = 24, outName = "supplemental_data_9", returnViln = F, 
          outDir = "./output/supplementalData/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
          assay = "RNA", min.pct = 0.25, only.pos = T)


seu.obj$celltype.l1 <- seu.obj$majorID_pertyName
seu.obj$celltype.l2 <- seu.obj$finalClusters
### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = "./output/supplementalData/supplemental_data_9_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "clusterID",
                                                   "Phase", "celltype.l1", "celltype.l2",  "name", "cellSource"), 
               skipEXPR = F, test = F,
               feats = c(
                   "PTPRC","CD3E","AIF1","CD68",
                   "CD34","ESAM","NOTCH3","ACTA2",
                   "FBLN1","COL1A1","IGFBP2","CXCL14",
                   "IL13RA2","DEFB1","DNASE1L3","CD55"
               )
              )

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
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.6, smallAxes = T) & NoLegend()
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



#column_title_gp = gpar(fontsize = 14, fontface = "bold")
# Export supplemental data
table(seu.obj$finalClusters, seu.obj$name) %>%
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
    write.csv(., "./output/supplementalData/syn_permRes.csv", row.names = F)

### DGE analysis using leinent Wilcoxon rank-sum test


### Fig 1d/e - DEG heatmap and scatter plot of DEGs in each major subset
#complete dge analysis within each major subset
seu.obj$finalClusters <- factor(gsub("/", "_", seu.obj$finalClusters))
seu.obj$finalClusters <- factor(seu.obj$finalClusters, levels = c("Synoviocyte_c0", "Synoviocyte_c1", "Synoviocyte_c2",
                                                                   "Macrophage_DC", "Synoviocyte_c4", "Synoviocyte_c5", 
                                                                   "Synoviocyte_c6", "Endothelial", "Fibroblast",
                                                                   "T cell", "Myofibroblast")) 

linDEG(seu.obj = seu.obj, groupBy = "finalClusters", comparision = "cellSource", outDir = paste0("./output/", outName,"/linDEG/"), outName = "allCells", cluster = NULL, labCutoff = 10, noTitle = F, labsHide = "^ENSECAG", contrast = c("OA", "Normal"),
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, returnPlots = F
                  )


#load dge results and create heatmap of degs
files <- list.files(path = "./output/allCells_syn/linDEG/", pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv, header = T)

cnts_mat <- do.call(rbind, df.list)  %>% 
    mutate(
        direction = case_when(
            avg_log2FC <= -1 ~ "Down",
            avg_log2FC >= 1 ~ "Up",
            avg_log2FC < 1 & avg_log2FC > -1 ~ "n.s."
        )
    ) %>% 
    filter(!direction == "n.s.") %>%
    group_by(cellType, direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = cellType, values_from = nRow) %>% 
    as.matrix() %>% t()

# Generate table with cell numbers for reviwer
as.data.frame(table(seu.obj$finalClusters, seu.obj$name)) %>% 
    pivot_wider(names_from = "Var2", values_from = "Freq") %>%
    rename(`Cell type` = Var1) %>%
    write.csv(., "./output/supplementalData/cell_counts_syn.csv", row.names = F)

colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"

orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

png(file = paste0("./output/", outName, "/",outName, "_deg_heat.png"), width=1500, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col = viridis(100), #circlize::colorRamp2(c(0,max(cnts_mat)), colors = c("white","red")),
              cluster_columns = F,
              column_title = gt_render(
                  paste0("<span style='font-size:18pt; color:black'># of DEGs</span><br>",
                         "<span style='font-size:12pt; color:black'>(OA vs Normal)</span>")
              ),
              show_column_names = TRUE,
              column_title_side = "top",
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
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()


### Supp fig xx -- heatmap of degs by each cluster
seu.obj$type <- factor(paste0(as.character(seu.obj$finalClusters), "--", as.character(seu.obj$cellSource)),
                       levels = paste0(rep(levels(seu.obj$finalClusters), each = 2), "--", c("Normal", "OA")))

files <- paste0("./output/allCells_syn/linDEG/", "allCells_", gsub(" ", "_", levels(seu.obj$finalClusters)), "_geneList.csv")
df.list <- lapply(files, read.csv, header = T)
res.df <- do.call(rbind, df.list) %>% filter(abs(avg_log2FC) > 1) %>% filter(!grepl("^ENS", X))

sig.mat <- matrix(nrow = length(unique(res.df$X)), ncol = length(levels(seu.obj$type)),
                  dimnames = list(unique(res.df$X),
                                  toupper(levels(seu.obj$type))))

for(i in 1:nrow(sig.mat)){
    for(j in 1:ncol(sig.mat)){
        cellType <- strsplit(colnames(sig.mat)[j], "--")[[1]][1]
        condition <- strsplit(colnames(sig.mat)[j], "--")[[1]][2]
        if(cellType %in% toupper(res.df[res.df$X == rownames(sig.mat)[i], ]$cellType)){
            lfc <- res.df[res.df$X == rownames(sig.mat)[i] & toupper(res.df$cellType) == cellType, ]$avg_log2FC
            if(lfc > 1 & condition == "OA"){
                sig.mat[i, j] <- "*"
            } else if(lfc < -1 & condition == "NORMAL"){
                sig.mat[i, j] <- "*"
            } else{
                sig.mat[i, j] <- ""
            }
        } else{
            sig.mat[i, j] <- ""
        }
    }
}

res.df <- res.df[!duplicated(res.df$X), ]

#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@counts)) #use raw count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

#filter matrix for DEGs and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% res.df$X]
mat_scaled <- t(apply(t(log1p(clusAvg_expression)), 1, scale))
colnames(mat_scaled) <- rownames(clusAvg_expression)
mat_scaled <- mat_scaled[ ,match(colnames(sig.mat), toupper(colnames(mat_scaled)))]
mat_scaled <- mat_scaled[match(rownames(sig.mat), rownames(mat_scaled)), ]  

#set annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name)
clus <- colz.base
names(clus) <- levels(seu.obj$finalClusters)
cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Normal","OA")

# heat_col <- viridis(option = "magma",100)
ha <- HeatmapAnnotation(
    Cluster = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
    Condition = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][2]})),
    border = TRUE,
    col = list(Cluster = clus, Condition = cond_colz)
)

#plot the data
ht <- Heatmap(
    mat_scaled,
    name = "mat",
    cluster_rows = F,
    row_title_gp = gpar(fontsize = 24),
    show_row_names = T,
    cluster_columns = F,
    top_annotation = ha,
    show_column_names = F,
    column_split = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
                          levels = levels(seu.obj$finalClusters)),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
            title = "Scaled expression",
            direction = "horizontal"
        ),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
    }
)

png(file = paste0("./output/", outName, "/", outName, "_fig3e.png"), width=3750, height=3750, res=400)
par(mfcol=c(1,1))   
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))

for(i in 1:length(levels(seu.obj$finalClusters))){
    decorate_annotation("Cluster", slice = i, {
        grid.text(paste0("c", (1:11) - 1)[i], just = "center")
    })
}
dev.off()


##################################
### END ALL CELL DATA ANALYSIS ###
##################################
