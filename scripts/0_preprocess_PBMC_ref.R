#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

#pbmc data

load10x(din = "./pbmc_data/", dout = "./output/s1_pbmc/", outName = "pbmc_data", testQC = F,
       nFeature_RNA_high = 3750, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 25000, nCount_RNA_low = 100, readCnts = T)

seu.obj <- sctIntegrate(din = "./output/s1_pbmc/", dout = "./output/s2/", outName = "2023_09_11_pbmc_data", vars.to.regress = "percent.mt", nfeatures = 2500)

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "2023_09_11_pbmc_data", test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "2023_09_11_pbmc_data", final.dims = 50, final.res = 0.4, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("CD3G", "TRAT1", "PRF1", "CTSW", 
                                     "MPEG1", "DRA", "CST3", "CD14", 
                                     "MS4A1", "CD79A", "LTC4S","TOP2A")
                       )

outName <- "allCells"

seu.obj <- readRDS("./output/s3/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds")
outName <- "pbmc"


namez <- lapply(colnames(seu.obj), function(x){
    strsplit(x, split = "\\.")[[1]][1]
})

unique(unlist(namez))


cellID.df <- read.csv("./eq_idents.csv")
# cellID$X <- NULL

# cellID.df <- cellID %>% mutate(eqName = case_when(sampleID == "s1" ~ "BEN",
#                                      sampleID == "s2" ~ "BOUNCE",
#                                     sampleID == "s3" ~ "IDO",
#                                     sampleID == "s4" ~ "JJ",
#                                     sampleID == "s5" ~ "MEMPHIS",
#                                     sampleID == "s6" ~ "ROANIE",
#                                     sampleID == "s7" ~ "STARDUST"),
#                   barcode = paste0(eqName,".",barcode)
#                  )

ct.list <- as.factor(cellID.df$clusterID_patel)
names(ct.list) <- cellID.df$eqName.barcode

seu.obj <- AddMetaData(seu.obj, metadata = ct.list, col.name = "clusterID_patel")


### Fig 1a: plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_patel",
           #    cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./patel_idents.csv", groupBy = "clusterID_patel", metaAdd = "celltype.l1")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./patel_idents.csv", groupBy = "clusterID_patel", metaAdd = "celltype.l2")


# meta.df <- seu.obj.sub@meta.data
# meta.df <- meta.df[ ,c(1,33,34,35)] %>% rownames_to_column(., "eqName.barcode") %>% mutate(eqName = strsplit(eqName.barcode, split = "\\.")[[1]][1],
#                                                                          sampleID = strsplit(orig.ident, split = "_")[[1]][2],
#                                                                          barcode = strsplit(eqName.barcode, split = "\\.")[[1]][2],
#                                                                                            sampleID_barcode = paste0(sampleID,"_",barcode)
#                                                                         )
# write.csv(meta.df[ ,c(8,1,9,2,7,6,3,4,5)],"./eq_idents.csv", row.names = F)

### Fig 1a: plot inital cluster umap
pi <- DimPlot(seu.obj.sub, 
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



### Fig 1a: plot inital cluster umap
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "celltype.l2",
           #    cols = colArray$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = F
 )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_UMAP_ctl2.png", sep = ""), width = 7, height = 7)

# saveRDS(seu.obj, file = "./output/s3/2023_09_11_pbmc_data_res0.4_dims50_dist0.5_neigh50_S3.rds")

seu.obj.tcell <- readRDS("/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/s3/Aug_10_2023_tcell_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")


seu.obj.sf <- readRDS("./output/s3/Sept_4_2023_allCells_2500_res0.9_dims45_dist0.3_neigh30_S3.rds")
seu.obj <- subset(seu.obj, invert = T,
                  subset = clusterID_patel == "NA")


