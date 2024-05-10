#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(scProportionTest)

########################################
### BEGIN MYELOID DATA PREPROCESSING ###
########################################

#load preprocessed allCells object
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./sf_idents_08-01-2023.csv", groupBy = "clusterID", metaAdd = "majorID")
outName <- "myeloid"
datE <- "Aug_8_2023"
nfeatures <- 2500

#subset on myeloid cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "myeloid")

#ensure subset was properly done
table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , preSub = T, nfeatures = nfeatures,
                      vars.to.regress = "percent.mt"
                       )

#check ideal clustering resolution
# seu.obj <- readRDS("./output/s2/Aug_8_2023_myeloid_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#dim reduction and unsupervised clustering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


### During analysis several low quality clusters were identified and need to removed 
### NOTE: a very small a cluster of fibroblasts was identified and is removed due to not belonging in this subset

#load in previously processed data and load metadata
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

#check cluster resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), test_dims = "40", algorithm = 3, prefix = "integrated_snn_res.")


#complete unsupervised clustering and dim reduction
seu.obj$clusterID_sub_orig <- seu.obj$clusterID_sub
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

### Yet another low quality subset was identified
### Remove and only repeat dim reduction due to how small the cluster is

#stash previous cluster IDs
seu.obj$clusterID_sub_orig2 <- seu.obj$clusterID_sub

#subset on inverted unwanted cells
seu.obj <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "12")

#complete unsupervised clustering and dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


######################################
### END MYELOID DATA PREPROCESSING ###
######################################

#Note: the processed Seurat object is now stored as a .rds file, "./output/s3/Sep_1_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds"

###################################
### BEGIN MYELOID DATA ANALYSIS ###
###################################

#load in processed data and metadata
seu.obj <- readRDS("./output/s3/Sep_1_2023_myeloid_2500_res0.4_dims40_dist0.3_neigh30_S3.rds")
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
seu.obj$majorID_sub_inc <- Idents(seu.obj)
seu.obj$majorID_sub_inc <- factor(gsub("_Macrophage", "_Mac", as.character(seu.obj$majorID_sub_inc)),
                                    levels = gsub("_Macrophage", "_Mac", levels(seu.obj$majorID_sub_inc)))
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub_inc, levels = levels(seu.obj$majorID_sub_inc)[c(8,9,7,5,4,3,2,6,1,10,12,11)])

#set higher level major idents
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "DC" , "1" = "Macrophage", 
                                   "2" = "Macrophage", "3" = "Macrophage",
                                   "4" = "Macrophage", "5" = "DC",
                                   "6" = "Monocyte", "7" = "Neutrophil",
                                   "8" = "Neutrophil", "9" = "DC",
                                   "10" = "DC", "11" = "DC"
                       ))
seu.obj$majorID_sub2 <- Idents(seu.obj)


#summary data for cell type percentages in myeloid cells
table(seu.obj$majorID_sub_inc, seu.obj$name) %>%
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
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = "supplemental_data_4", returnViln = F, 
                      outDir = paste0("./output/supplementalData/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T)


### Data supplemental - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, dir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",datE,outName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID", "clusterID_sub", "majorID_sub", "name", "cellSource"), 
               skipEXPR = F, test = F,
                           feats = c("CD68", "FLT3", "CD1C")
                          
                          )

#load in cell type colors
colArray <- read.csv("./sf_idents_08-01-2023.csv")
colArray.sub <- colArray[colArray$majorID == "myeloid",]

### Fig 2a - Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              cols = colArray.sub$newCol,
              label = T,
              label.box = T,
              shuffle = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.6, smallAxes = TRUE) 
ggsave(paste0("./output/", outName, "/", outName, "_fig2a.png"), width = 7, height = 7)


### Fig extra - Create UMAP by majorID
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
ggsave(paste("./output/", outName, "/", outName, "_majorSub_UMAP.png", sep = ""), width = 7, height = 7)


### Fig extra - plot key feats
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


### Fig 2b - Create violin plots for key feats
features <- c("AIF1","DRA","ITGAX",
              "SELP", "CSF3R","SELL",
              "CXCL8","VCAN", "LYZ",
              "AQP9",
              "MSR1","CSF1R",
              "CCL2","CCL8","LYVE1",
              "TPPP3",
              "MARCO","CD5L","GPNMB",
              "FLT3",
              "BATF3","DNASE1L3","IRF8",
              "CD1C","FCER1A",
              "TCF4","MS4A1",
              "LAMP3","CCR7",
              "TOP2A")

pi <- VlnPlot(object = seu.obj,
              pt.size = 0,
              group.by = "majorID_sub",
              combine = T,
              stack = T,
              cols = colArray.sub$newCol[c(8,9,7,5,4,3,2,6,1,10,12,11)],
              fill.by = "ident",
              flip = T,
              features = rev(features)
             ) & NoLegend() & theme(axis.ticks.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(size = 7),
                                    axis.text = element_text(size = 7),
                                    axis.line = element_blank(),
                                    strip.text.y.right = element_text(size = 7, hjust = 0)
                                   ) + theme(plot.margin = margin(c(3,0,0,10)))
ggsave(paste("./output/", outName, "/", outName, "_fig2b.png", sep = ""), width = 4, height = 6)


### Fig extra: make dot plots of key features
pi <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
               features = features
              ) + theme(axis.title = element_blank(),
                        axis.text = element_text(size = 12)
                       )
ggsave(paste("./output/", outName, "/", outName, "_majorDot.png", sep = ""), width = 10, height = 5)


### Fig supp 2a - M1 vs M2 signatures
#set and calc modules
modulez <- list("Pro-inflammatory" = c("AZIN1", "CD38","CD86","CXCL10","FPR2","GPR18","IL12B","IL18","IRF5","NFKBIZ","NOS2","PTGS2","TLR4","TNF"),
                "Anti-inflammatory" = c("ALOX15", "ARG1", "CHIL3", "CHIL4","EGR2", "IL10","IRF4","KLF4","MRC1","MYC","SOCS2","TGM2")
                )
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
features <- names(modulez)

modulez <- c(list("Enrichment score" = names(modulez)), modulez)
plots <- lapply(modulez, function(x){
    majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
             features = rev(unname(unlist(x)))
            ) + coord_flip() + theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     legend.position = "right",
                                     legend.direction = "vertical",
                                     axis.title = element_blank(),
                                     plot.margin = margin(3, 0, 3, 0, "pt")
                                    ) + scale_colour_viridis(option="magma", name='Average\nexpression', limits = c(-2.5,2.5)) + guides(size = guide_legend(nrow = 3))
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

plots$`Enrichment score` <- plots$`Enrichment score` + theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0)) + scale_y_discrete(position = "right") + scale_colour_viridis() + guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))
p <- Reduce( `+`, plots ) +  plot_layout(guides = "collect", design = patch, 
                                         height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length)))))
ggsave(paste("./output/", outName, "/", outName, "_supp2a.png", sep = ""), width = 6, height = 7)


### Fig supp 2b umap by sample
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
p <- formatUMAP(pi) + NoLegend() + theme(legend.position = "top", 
                                         legend.justification = "center",
                                         legend.direction = "horizontal",
                                         legend.title = element_text(size = 12)
                                        ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_supp2b.png", sep = ""), width = 7, height = 7)


### Fig extra - stacked bar graph by majorID_sub
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "majorID_sub") +
scale_fill_manual(labels = levels(seu.obj$name), values = levels(seu.obj$colz)) + 
theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12)) 
ggsave(paste("./output/", outName,"/",outName, "_stackedBar.png", sep = ""), width = 8, height = 5)


### Fig extra - evaluate cell frequency by majorID_sub using wilcoxon rank sum test
freqy <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name", showPval = T,
               namez = "name", 
               colz = "colz"
              )
ggsave(paste("./output/", outName, "/",outName, "_freqPlots_majorID_sub.png", sep = ""), width = 8.5, height = 9)


### Fig 2c: monte-carlo permitation test
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
ht <- daOR(seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "name", 
           outName = outName, outDir = paste0("./output/", outName, "/"), 
           t_map = TRUE, cluster_order = p$data$clusters)

png(file = paste0("./output/", outName, "/", outName, "_OR_heat.png"), width = 2000, height = 1000, res = 400)
par(mfcol = c(1, 1))         
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "left")
dev.off()

# Export supplemental data
table(seu.obj$majorID_sub, seu.obj$name) %>%
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
    write.csv(., "./output/supplementalData/myeloid_permRes.csv", row.names = F)

### Extra analysis - evaluate differntial anundance using neighborhood based miloR analysis
library(miloR)
library(BiocParallel)

runMilo <- function(
    seu.obj = NULL, 
    da_design = NULL, 
    subName = "", 
    blocked = TRUE, 
    ...
    ){
    
    # Hard coded metadata re-naming
    seu.obj$Sample <- seu.obj$name
    seu.obj$Condition <- seu.obj$cellSource

    # Convert from Seurat to sce object
    sce <- as.SingleCellExperiment(seu.obj)
    reducedDim(sce, "PCA") <- seu.obj@reductions$pca@cell.embeddings
    reducedDim(sce, "UMAP") <- seu.obj@reductions$umap@cell.embeddings

    # Preprocess using miloR to ID neighboorhoods
    milo.obj <- Milo(sce)
    milo.obj$Sample <- droplevels(factor(milo.obj$Sample))
    milo.obj <- buildGraph(milo.obj, k = 30, d = 40)
    milo.obj <- makeNhoods(milo.obj, prop = 0.2, k = 30, d = 40, refined = TRUE, refinement_scheme = "graph")
    p <- plotNhoodSizeHist(milo.obj)
    ggsave(paste0("./output/", outName, "/", outName, "_NhoodSize.png"), width = 7, height = 7)
    
    milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), samples = "Sample")

    # Set up metadata
    rownames(da_design) <- da_design$Sample
    da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]

    # Calc distance between neighborhoods and test for DA
    milo.obj <- calcNhoodDistance(milo.obj, d = 40)
    if(blocked){
        da_results <- testNhoods(milo.obj, design = ~ Batch + Condition, design.df = da_design)
    } else{
        da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
    }
    
    n_diff <- da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.1) %>% nrow()
    if(n_diff == 0){
        message(
            paste(
                "No differentially abundant Nhoods at alpha = 0.1!",
                "The lowest spaitally adjusted P value is:", min(da_results$SpatialFDR),
                "\n Although not reccomended, you can increase alpha to the lowest",
                "spaitally adjusted P value to get an idea of which regoins of the",
                "UMAP are trendy."
            )
        )
    }

    # Plot the results (by neighborhood)
    milo.obj <- buildNhoodGraph(milo.obj)
    p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods=!is.na(da_results$logFC), ...)
    ggsave(paste("./output/", outName, "/", subName, "_milo.png", sep = ""), width = 6, height = 6)
    return(list(p, milo.obj))
}

# Set up metadata
da_design <- as.data.frame(list(
    "Sample" = factor(c("Normal_1", "Normal_2", "Normal_3", "OA_1", "OA_2", "OA_3")),
    "Condition" = factor(c("Normal", "Normal", "Normal", "OA", "OA", "OA")),
    "Batch" = factor(c(1, 2, 3, 1, 2, 4))
))

p <- runMilo(seu.obj = seu.obj, da_design = da_design, subName = "OA_vs_Norm", blocked = T, alpha = 0.1)

p0 <- p[[1]] + ggtitle("OA versus Normal") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste0("./output/", outName, "/", outName, "_milo_sig.png"), width = 7, height = 7)

milo.obj <- p[[2]] 
rownames(da_design) <- da_design$Sample
da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
da_results <- testNhoods(milo.obj, design = ~ Batch + Condition, design.df = da_design)
da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.95) %>% nrow()

p0 <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods = !is.na(da_results$logFC), alpha = 1)
p0 <- p0 + ggtitle("OA versus Normal") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste("./output/", outName, "/", outName, "_milo_all.png", sep = ""), width = 7, height = 7)


### Fig supp 2c - genes that define skewed DCs
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "majorID_sub", idents.1 = "cDC2", idents.2 = "cDC1", bioRep = "name", padj_cutoff = 0.01, lfcCut = 1, strict_lfc = T,
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "cDC2_VS_cDC1", idents.1_NAME = "cDC2", idents.2_NAME = "cDC1", returnVolc = T, doLinDEG = F, paired = T, lowFilter = T, dwnSam = F, setSeed = 24, topn = c(15,15), addLabs = ""
                    )
p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c0 (cDC2)", leftLab = "Up in c5 (cDC1)", lfcCut = 1) + labs(x = "log2(FC) cDC2 vs cDC1") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_c0vc5_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig 2d: DGE using wilcoxon
seu.obj$allCells <- "DGE analysis of myeloid cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("OA", "Normal"))
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/fig2d_"), 
       outName = outName, labCutoff = 10, saveGeneList = T,
       logfc.threshold = 0.58, pValCutoff = 0.01
                  )


### Fig 2d - complete pseudobulk DGE by all cells
seu.obj$allCells <- "All_myeloid_cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
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
          paired = T, pairBy = "horse", strict_lfc = F, labSize = 4.5,
          idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

pi  <- prettyVolc(plot = p[[1]], rightLab = "Up in OA", leftLab = "Up in Normal", arrowz = T, lfcCut = 1) + labs(x = "log2(Fold change)") + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank(),
                                                                                                                                  plot.title = element_text(face = "bold", hjust = 0.5, size = 20, vjust = 2)) + ggtitle("Myeloid cells (OA vs Normal)")
ggsave(paste("./output/", outName, "/", outName, "_fig1e.png", sep = ""), width = 7, height = 7)


#confirm overlap of key DEGS using both approaches
lindeg.df <- read.csv("./output/myeloid/fig2d_myeloid_DGE_analysis_of_myeloid_cells_geneList.csv")
pb.df <- read.csv("./output/myeloid/pseudoBulk/allCells/myeloid_cluster_allCells_all_genes.csv") %>% filter(abs(log2FoldChange) > 1)
table(pb.df$gene %in% lindeg.df$X)
# FALSE  TRUE 
#     6     7 

### By cluster
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("./output/", outName, "/pseudoBulk/"), min.cell = 5,
         grepTerm = "Normal", grepLabel = c("Normal", "OA")
)

df <- read.csv(paste0("./output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"), row.names = 1)
df$horse <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
df <- df %>% mutate(horse = ifelse(sampleID == "Normal_3", "horse_4", paste0("horse_",horse)))
write.csv(df, paste0("./output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("./output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
          padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName, "/pseudoBulk/"), 
          outName = outName, returnVolc = F, 
          paired = T, pairBy = "horse", strict_lfc = F,
          idents.1_NAME = "OA", idents.2_NAME = "Normal",
          inDir = paste0("./output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

# Export supp data
files <- lapply(c("All_myeloid_cells", levels(seu.obj$majorID_sub)), function(x){
    list.files(path = paste0("./output/myeloid/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

files <- unlist(files)
df.list <- lapply(files, read.csv, header = T)
supp.df <- do.call(rbind, df.list) %>% filter(abs(log2FoldChange) > 1)
write.csv(supp.df, "./output/supplementalData/supplemental_data_6_gene_list.csv", row.names = F)


### Supp fig xx - heatmap of dge results by major cell types -- number of DEGS
files <- lapply(levels(seu.obj$majorID_sub), function(x){
    list.files(path = paste0("./output/myeloid/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

files <- unlist(files)
df.list <- lapply(files, read.csv, header = T)

cnts_mat <- do.call(rbind, df.list)  %>% 
    mutate(
        direction = case_when(
            log2FoldChange >= 1 ~ "Up",
            log2FoldChange <= -1 ~ "Down",
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


cnts_mat[is.na(cnts_mat)] <- 0
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

png(file = paste0("./output/", outName, "/",outName, "_fig1f.png"), width=1500, height=2000, res=400)
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
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 14, col = "black"),
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                  if(cnts_mat[i, j] > 5) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
                  } else if(cnts_mat[i, j] <= 5) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "grey80"))
                  }                          
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"), show_heatmap_legend = FALSE)
dev.off()


### Supp fig xx -- heatmap of degs by each cluster
# levels(seu.obj$majorID_sub_inc) <- paste0("(c", ((1:12) - 1), ") ", levels(seu.obj$majorID_sub_inc))
seu.obj$type <- factor(paste0(as.character(seu.obj$majorID_sub_inc), "--", as.character(seu.obj$cellSource)),
                       levels = paste0(rep(levels(seu.obj$majorID_sub_inc), each = 2), "--", c("Normal", "OA")))


res.df <- do.call(rbind, df.list) %>% filter(abs(log2FoldChange) > 1)

sig.mat <- matrix(nrow = length(unique(res.df$gene)), ncol = length(levels(seu.obj$type)),
                  dimnames = list(unique(res.df$gene),
                                  toupper(levels(seu.obj$type))))
# colnames(sig.mat) <- gsub("(.*.) ", "", colnames(sig.mat))

for(i in 1:nrow(sig.mat)){
    for(j in 1:ncol(sig.mat)){
        cellType <- strsplit(colnames(sig.mat)[j], "--")[[1]][1]
        condition <- strsplit(colnames(sig.mat)[j], "--")[[1]][2]
        if(cellType %in% res.df[res.df$gene == rownames(sig.mat)[i], ]$gs_base){
            lfc <- res.df[res.df$gene == rownames(sig.mat)[i] & res.df$gs_base == cellType, ]$log2FoldChange
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

res.df <- res.df[!duplicated(res.df$gene), ]


#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@counts)) #use raw count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

#filter matrix for DEGs and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% res.df$gene]
mat_scaled <- t(apply(t(log1p(clusAvg_expression)), 1, scale))
colnames(mat_scaled) <- rownames(clusAvg_expression)
mat_scaled <- mat_scaled[ ,match(colnames(sig.mat), toupper(colnames(mat_scaled)))]
mat_scaled <- mat_scaled[match(rownames(sig.mat), rownames(mat_scaled)), ]  

#set annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name2)
clus <- colArray.sub$newCol
names(clus) <- paste0("(c", ((1:12) - 1), ") ", levels(seu.obj$majorID_sub_inc))
cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Normal","OA")

ha <- HeatmapAnnotation(
    Cluster = factor(paste0("(c", rep(((1:12) - 1), each = 2), ") ", unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]}))),
                     levels = paste0("(c", ((1:12) - 1), ") ", levels(seu.obj$majorID_sub_inc))),
    Condition = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][2]})),
    border = TRUE,
    col = list(Cluster = clus, Condition = cond_colz),
    annotation_legend_param = list(
        Cluster = list(direction = "horizontal", nrow = 3),
        Condition = list(direction = "vertical", nrow = 2)
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

lgd1 <- Legend(labels = paste0("(c", ((1:12) - 1), ") ", levels(seu.obj$majorID_sub_inc)),
               legend_gp = gpar(fill = clus), 
               title = "Clusters", 
               direction = "vertical",
               nrow = 3, 
               gap = unit(0.6, "cm")
              )

lgd2 <- Legend(labels = names(cond_colz),
               legend_gp = gpar(fill = cond_colz), 
               title = "Condition", 
               direction = "vertical",
               nrow = 2
              )

pd <- packLegend(lgd1, lgd2, max_width = unit(45, "cm"), 
    direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(0.5, "cm"))

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
                          levels = levels(seu.obj$majorID_sub_inc)),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
        title = "Average expression",
        direction = "horizontal"
        ),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
    }
)

png(file = paste0("./output/", outName, "/", outName, "_fig3e.png"), width=3500, height=3500, res=400)
par(mfcol=c(1,1))   
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom", 
     annotation_legend_list = pd, annotation_legend_side = "top")

for(i in 1:length(levels(seu.obj$majorID_sub_inc))){
    decorate_annotation("Cluster", slice = i, {
        grid.text(paste0("c", (1:12) - 1)[i], just = "center")
    })
}
dev.off()


### Fig supp 2d - run/plot gsea results
dgea.df <- read.csv("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/myeloid/fig2d_myeloid_DGE_analysis_of_myeloid_cells_geneList.csv")
geneListUp <- dgea.df %>% arrange(p_val_adj) %>% filter(p_val_adj < 0.01, avg_log2FC > 0) %>% pull(X)
geneListDwn <- dgea.df %>% arrange(p_val_adj) %>% filter(p_val_adj < 0.01, avg_log2FC < 0) %>% pull(X)

p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5", upCol = "red", dwnCol = "blue", upOnly = T, termsTOplot=35, trimTerm = T)+ theme(axis.title = element_text(size = 16))
pi <- p + scale_x_continuous(limits = c(-50,ceiling(max(p$data$x_axis)*1.05)), 
                             breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
                             name = "-log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("./output/", outName, "/", outName, "_supp2d.png"), width = 7.5, height = 9)


### Fig supp 2e - use enrichment scoring to identify which cells are enriched in the GSEA terms
enriched <- pi$data
enriched <- enriched %>% mutate(labz = ifelse(grepl("MIGRATION|CHEMOTAXIS|ACTIVATION",Description),Description,NA))
enriched$geneID <- as.list(strsplit(enriched$geneID, "/"))
modulez <- enriched[!is.na(enriched$labz),]$geneID 

#organize the data
names(modulez) <- enriched[!is.na(enriched$labz),]$Description

#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched')) + coord_flip()
ggsave(paste("./output/", outName, "/", outName, "_supp2e.png", sep = ""), width = 9, height = 5)


### Fig supp 2e - plot enrichment scores for gsea terms by cell type
enriched$geneID <- as.list(strsplit(enriched$geneID, "/"))
modulez <- enriched[!is.na(enriched$labz),]$geneID 

#organize the data
names(modulez) <- enriched[!is.na(enriched$labz),]$Description

#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_",term,"_dots_celltypes.png", sep = ""),width = 10,height=10)


### Fig 2e - split feature plots highlighting key degs
set.seed(12)
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Normal", "OA"))
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
# features <- c("SPP1","","CD1A3","IL1RN","ACTA2","BCAM","KLRK1", "SYTL3") #supplemtnal features
features <- c("GPNMB","CPM","SPP1","TIMP1","TREM1","MARCO","SAMD9L")
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


ggsave(paste("./output/", outName, "/", outName, "_splitFeats.png", sep = ""), width = 10, height = 4)


legg <- FeaturePlot(seu.obj.sub,features = features[1]) + 
    theme(
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center",
    ) + scale_colour_viridis(breaks = c(0, 1), 
                             limits = c(0, 1),
                             label = c("low", "high"),
                             option = "magma")
legg <- get_legend(legg)
ggsave(plot = legg, paste("./output/", outName, "/", outName, "_splitFeats.png", sep = ""), width = 10, height = 4)

#################################
### END MYELOID DATA ANALYSIS ###
#################################
