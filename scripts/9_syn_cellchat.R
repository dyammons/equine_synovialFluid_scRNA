#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(CellChat)
outName <- "cellchat"
subName <- "syn"

### Load final data
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

seu.obj.backUp <- seu.obj


### run cell chat on treated
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "OA")

cnts <- seu.obj@assays$RNA@data
#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                     gene_input = "rownames", 
                                     gene_output = "rownames", 
                                     input_species = "horse",
                                     output_species = "human",
                                     non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
meta$chatID <- meta$finalClusters
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "chatID") 

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/cellChatobj_oa_finalClusters_231228.rds")
# saveRDS(cellchat, "./output/cellchat/cellChatobj_oa_clusID.rds")


#run cell chat on transfered naive
seu.obj <- seu.obj.backUp

Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Normal")

table(seu.obj$cellSource)

cnts <- seu.obj@assays$RNA@data
#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                     gene_input = "rownames", 
                                     gene_output = "rownames", 
                                     input_species = "horse",
                                     output_species = "human",
                                     non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
meta$chatID <- meta$finalClusters
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "chatID") 

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/cellChatobj_norm_finalClusters_231228.rds")
# saveRDS(cellchat, "./output/cellchat/cellChatobj_norm_clusID.rds")



#### do the comparision analysis  --- this analysis may be confounded by differences in filter btwn tx and naive samples
outName <- "cellchat"
subName <- "syn"

cellchat.norm <- readRDS("./output/cellchat/cellChatobj_norm_finalClusters_231228.rds")
cellchat.oa <- readRDS("./output/cellchat/cellChatobj_oa_finalClusters_231228.rds")

#load in the cell type colors
colz.base <- colz.base <- c("#00C1A7","#00A6FF","#00BADE", "#64B200", "#00B5ED", "#00C0BB", "#619CFF", "#AEA200", "#DB8E00", "#B385FF", "#F8766D") 
names(colz.base) <- unique(seu.obj$finalClusters)[c(1,3,2,5,7,4,6,10,8,9,11)]

# cellchat.norm <- readRDS("./output/cellchat/cellChatobj_norm_clusID.rds")
# cellchat.oa <- readRDS("./output/cellchat/cellChatobj_oa_clusID.rds")

object.list <- list(Normal = cellchat.norm, OA = cellchat.oa)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP") 
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_int.png")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + 
    theme(title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 16)
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + coord_flip()
ggsave(outfile, height=2,width=6)


#create heatmap using CellChat built-in function
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         cluster.rows = F,
                         cluster.cols = F, color.use = colz.base)

#extract the data and customize the heatmap


colz.base <- colz.base[match(rownames(gg2@matrix), names(colz.base))]
ha <- HeatmapAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )
row_ha <- rowAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )

outfile <- paste0("./output/", outName, "/", subName, "_fig5d.png")
png(file = outfile, width=2250, height=2500, res=500)

cusHeat <- Heatmap(gg2@matrix, na_col = "white", col = gg2@matrix_color_mapping, 
                   bottom_annotation = ha, left_annotation = row_ha,
                   show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                   column_title = "Differenital interaction strength", column_title_gp = gpar(fontsize = 10), column_names_rot = 60,
                   row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                   cluster_rows=T, cluster_columns=T,
                  
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                               border = NA,
                                               legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
                 )

cusHeat
dev.off()



#plot information flow
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_infoFLOW.png")
png(file = outfile, width=1000, height=500)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = colz.base)
}

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_int2D.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,2), xpd=TRUE)
gg[[1]] + gg[[2]]
dev.off()




#extract data and customize plots
gg1.df <- gg[[1]]$data
gg1.df$data_type <- "Normal"

gg2.df <- gg[[2]]$data
gg2.df$data_type <- "OA"

gg.df <- rbind(gg1.df,gg2.df)

colz.df <- as.data.frame(colz.base)
colz.df$labels <- rownames(colz.df)

gg.df <- gg.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

#set cell types to label
gg.df <- gg.df %>% group_by(data_type) %>% arrange(desc(Count)) %>% mutate(lab=ifelse(row_number() <= 5, as.character(labels), NA)) %>% ungroup()

#make fig 4b
pis <- lapply(c("Normal","OA"),function(z){
    gg.df.sub <- gg.df %>% filter(data_type == z)

    ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, colour = labels, label=lab)) + 
            ggtitle(z) +
            geom_point() + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black") + 
            theme_classic() + 
    theme(axis.title = element_text(size= 10),
          axis.text = element_text(size= 8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
                 ) + NoLegend()
})

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 2, guides = 'collect') & scale_size_continuous(limits = c(min(gg.df$Count), max(gg.df$Count))) & xlim(0, 10) & ylim(0, 10) & scale_color_manual(values = unname(colz.base), labels = names(colz.base))
ggsave(pi2, file = paste0("./output/", outName, "/", outName, "_fig5e.png"), width = 6.25, height = 3)


### Fig supp 5 - unstacked bargraph of network information flow

#plot information flow using CellChat built-in function
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

#extract data and plot using custom code
res.df <- gg2$data %>% filter(pvalues < 0.05)

orderz <- res.df[res.df$group == "Normal", ] %>% left_join(res.df[res.df$group == "OA", ], by = "name") %>% mutate(pct = contribution.x/contribution.y) %>% pull(name) %>% as.character()

res.df$name <- factor(res.df$name, levels = orderz)
res.df$group <- factor(res.df$group, levels = c("Normal", "OA"))
p <- ggplot(res.df, aes(x = contribution.scaled, y = name, fill = factor(group))) +
            geom_bar(stat = "identity", position = "dodge", width = 1, colour="white") +
            theme_classic() +
            theme(title = element_text(size= 14),
                  legend.title = element_blank(),
                  legend.text = element_text(size= 12),
                  legend.position = "top",
                  legend.direction = "horizontal",
                  plot.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key.size = unit(1,"line"),
                  axis.title = element_text(size = 18),
                  axis.text = element_text(size = 16),
                  plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
            ) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0,0)) + 
            ylab(label = "Interaction network") +
            xlab(label = "Relative information flow") + 
            guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + NoLegend()

ggsave(p, file = paste0("./output/", outName, "/", subName, "_supp7a.png"))



print(cellchat.oa@netP$pathways[cellchat.oa@netP$pathways %in% cellchat.norm@netP$pathways])
print(cellchat.norm@netP$pathways[!cellchat.norm@netP$pathways %in% cellchat.oa@netP$pathways])



### Fig 8d/8e + Fig supp 8b-d: circos plots & violin plots

groupz.df <- table(seu.obj$majorID, seu.obj$finalClusters) %>% melt() %>% filter(value > 0)
groupzNames <- groupz.df$Var.1
names(groupzNames) <- groupz.df$Var.2


pathwayz <- c("CXCL")
pathwayz <- c("PDGF", "NOTCH")
pathways <- pathwayz
subName <- ""

cellchat <- cellchat.oa
cellchat <- cellchat.norm

colz.base <- colz.base <- c("#00C1A7","#00A6FF","#00BADE", "#64B200", "#00B5ED", "#00C0BB", "#619CFF", "#AEA200", "#DB8E00", "#B385FF", "#F8766D") 
names(colz.base) <- unique(seu.obj$finalClusters)[c(1,3,2,5,7,4,6,10,8,9,11)]

keyList <- lapply(pathwayz, function(pathway){
    
    #extract required plotting data
    lrData <- as.data.frame(cellchat@LR)
    net <- cellchat@net
    prob <- net$prob
    prob <- prob[,,rownames(lrData[lrData$LRsig.pathway_name == pathway,])]
    prob.sum <- apply(prob, c(1,2), sum)

    #identify which cell types are active in pathway
    idx1 <- which(Matrix::rowSums(prob.sum) == 0)
    idx2 <- which(Matrix::colSums(prob.sum) == 0)
    idx <- intersect(idx1, idx2)
    net <- prob.sum[-idx, ]
    if(is.matrix(net)){
        net <- net[, -idx]
        cellTypeOFinterest <- rownames(net)
    }else{
        net <- net[-idx]
        cellTypeOFinterest <- names(net)
    }
    
    #conditional to fix problem that occurs if all cell types are active in the plot
    if(is.null(cellTypeOFinterest)){
        cellTypeOFinterest <- levels(cellchat@idents)
    }

    #grey out cell types not involved
    colz2 <- colz.base
    colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("./output/", outName, "/", subName, "/", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate_mod(cellchat, layout = "chord", signaling = pathway, color.use = colz2, remove.isolate = F, big.gap = 5, group = groupzNames, title.name = NULL)
    dev.off()

    #get the active features in the pathway
    genez <- lapply(pathways, function(x){extractEnrichedLR(cellchat, signaling = x, geneLR.return = TRUE, enriched.only = T)[["geneLR"]]})
    names(genez) <- pathways
    
    keyList <- list("genez" = genez, "cellTypeOFinterest" = cellTypeOFinterest)
    return(keyList)
})

mapped_genes <- orthogene::map_genes(genes = "PDGFA",
                                     species = "horse", 
                                     drop_na = FALSE)

#manually plot
Idents(seu.obj) <- "finalClusters"
seu.obj.sub <- subset(seu.obj, idents = keyList[[1]]$cellTypeOFinterest)
seu.obj.sub$finalClusters <- as.factor(seu.obj.sub$finalClusters)

#plot expression using Seurat function
pi <- VlnPlot(
    object = seu.obj.sub,
    pt.size = 0,
    same.y.lims = T,
    flip = T,
    group.by = "finalClusters",
    fill.by = "ident",
    split.by = "cellSource",
#         cols = colz2,
    stack = T,
    combine = FALSE,
    features = c("ENSECAG00000038632", keyList[[1]]$genez$PDGF[2])
) + NoLegend() + theme(axis.title.x = element_blank(),
                       axis.title.y.right = element_blank()) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + ggtitle("PDGF network")

ggsave(paste("./output/", outName, "/", subName ,"supp7b.png", sep = ""), height = 3, width = 6)



#manually plot
Idents(seu.obj) <- "finalClusters"
seu.obj.sub <- subset(seu.obj, idents = keyList[[2]]$cellTypeOFinterest)
seu.obj.sub$finalClusters <- as.factor(seu.obj.sub$finalClusters)

#plot expression using Seurat function
pi <- VlnPlot(
    object = seu.obj.sub,
    pt.size = 0,
    same.y.lims = T,
    flip = T,
    group.by = "finalClusters",
    fill.by = "ident",
    split.by = "cellSource",
#         cols = colz2,
    stack = T,
    combine = FALSE,
    features = keyList[[1]]$genez$NOTCH
) + NoLegend() + theme(axis.title.x = element_blank(),
                       axis.title.y.right = element_blank()) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + ggtitle("NOTCH network")

ggsave(paste("./output/", outName, "/", subName ,"supp7c.png", sep = ""), height = 3, width = 6)








