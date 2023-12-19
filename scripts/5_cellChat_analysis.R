#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(CellChat)

######################################
### RUN CELLCHAT ON ANNOTATED DATA ###
######################################

#load in the annotated dataset
# seu.obj <- readRDS("./output/s3/Oct_30_2023_allCells_annotated.rds")
seu.obj <- readRDS("./output/s3/Nov_29_2023_allCells_annotated.rds")
seu.obj.backup <- seu.obj

### run cell chat on OA samples
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "OA")

cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "horse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
# meta$clusterID_sub <- paste("c_",meta$clusterID_sub, sep="")
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2") #  usin gpredicted ct.l2

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# saveRDS(cellchat, "./output/cellchat/Oct_30_2023_cellChatobj_oa_ctl2.rds")
saveRDS(cellchat, "./output/cellchat/Nov_30_2023_cellChatobj_oa_ctl2.rds")


#run cell chat on transfered naive

seu.obj <- seu.obj.backup

Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Normal")

cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "horse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
# meta$clusterID_sub <- paste("c_",meta$clusterID_sub, sep="")
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2") #  using predicted ct.l2

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# saveRDS(cellchat, "./output/cellchat/Oct_30_2023_cellChatobj_norm_ctl2.rds")
saveRDS(cellchat, "./output/cellchat/Nov_30_2023_cellChatobj_norm_ctl2.rds")



#### do the comparision analysis  --- this analysis may be confounded by differences in filter btwn tx and naive samples
outName <- "cellchat"
subName <- ""

seu.obj <- seu.obj.backup

# cellchat.norm <- readRDS("./output/cellchat/Oct_30_2023_cellChatobj_norm_ctl2.rds")
# cellchat.oa <- readRDS("./output/cellchat/Oct_30_2023_cellChatobj_oa_ctl2.rds")
cellchat.norm <- readRDS("./output/cellchat/Nov_30_2023_cellChatobj_norm_ctl2.rds")
cellchat.oa <- readRDS("./output/cellchat/Nov_30_2023_cellChatobj_oa_ctl2.rds")

#get overlapping pathways
cellchat.oa@netP$pathways[cellchat.oa@netP$pathways %in% cellchat.norm@netP$pathways]


object.list <- list(Normal = cellchat.norm, OA = cellchat.oa)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_int.png")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + 
    theme(title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          #axis.title.y = element_blank(),
          axis.text = element_text(size = 16)
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1"))  + scale_x_discrete(expand = c(0,0), breaks = c("Normal","OA")) + coord_flip()
#gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
#gg1 + gg2
ggsave(outfile, height=2,width=6)


data.df <- compareInteractions(cellchat, show.legend = F)$data
data.df

data.df$Dataset <- factor(data.df$dataset, levels = c("OA", "Normal"))
p <- ggplot(data.df, aes(x=count, y=Dataset, fill=dataset)) + geom_bar(stat="identity") + theme_void() + 
    theme(title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 16, hjust = 1),
          axis.text.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
          plot.margin = margin(4, 14, 4,4)
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1"))  + scale_y_discrete(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + NoLegend() + xlab("Number of inferred interactions")
ggsave(outfile, height=1,width=6)

levels(cellchat@idents) <- names(colz.base)
#pretty worthless
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_net.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = colz.base)
}

outfile <- paste0("./output/", outName,"//",outName,"_cellchat_cieVh_int2D.png")
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
            geom_point(color = gg.df.sub$colz.base) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black", box.padding = 1,  min.segment.length = 0.25) + 
            theme_classic() + 
    theme(axis.title.x = element_text(size= 10),
          axis.title.y = element_text(size= 10),
          axis.text = element_text(size=8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
                 ) + NoLegend()
})

pis[1][[1]] <- pis[1][[1]] + theme(axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()) 

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 1, guides = 'collect') & scale_size_continuous(limits = c(min(gg.df$Count), max(gg.df$Count))) & xlim(0, 12.5) & ylim(0, 12.5) & scale_color_manual(values = unname(colz.base), labels = names(colz.base))
# pi <- pi1 + pi2 + plot_layout(ncol = 2, widths = c(0.22,0.78)) 
ggsave(pi2, file = paste0("./output/", outName, "/", outName, "_interactionScater.png"), width = 3.25, height = 6)




#make heatmap -- pretty cool
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_weight.png")
png(file = outfile, width=2500, height=2500, res=500)
# gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T, color.use = colz.base)
gg2
dev.off()

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_weight.png")
png(file = outfile, width=2800, height=2500, res=500)


cusHeat <- Heatmap(gg2@matrix, na_col = "white", col = gg2@matrix_color_mapping, bottom_annotation = gg2@bottom_annotation, left_annotation = gg2@left_annotation, 
                  show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = "left",row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                  column_title = "Differenital interaction strength", column_title_gp = gpar(fontsize = 10), column_names_rot = 60,
                  row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, #at = colorbar.break,
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
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

res.df <- gg1$data %>% filter(pvalues < 0.05)

orderz <- res.df[res.df$group == "Normal", ] %>% left_join(res.df[res.df$group == "OA", ], by = "name") %>% mutate(pct = contribution.x/contribution.y) %>% pull(name) %>% as.character()

res.df$name <- factor(res.df$name, levels = orderz)
res.df$group <- factor(res.df$group, levels = c("Normal", "OA"))
p <- ggplot(res.df, aes(x = contribution, y = name, fill = factor(group))) +
    geom_bar(stat = "identity", position = "fill", width = 1, colour="white") +
    theme_classic() + geom_vline(xintercept = 0.5, linetype="dashed", color = "grey50", size=1) +
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

ggsave(p, file = paste0("./output/", outName, "/", subName, "/cellchat_oaVn_infoFLOW.png"))








###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}



outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_out.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,1), xpd=TRUE)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_in.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,1), xpd=TRUE, res = 500)
i <- 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_in.png")
png(file = outfile, width=1000, height=500)
pathways.show <- c("CD40") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  plot <- netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

plot
dev.off()



#make heatmap -- pretty cool
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_cnt.png")
png(file = outfile, width=2500, height=2500, res=500)
# gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "count",
                         cluster.rows = T,
                         cluster.cols = T)
gg2
dev.off()

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_cnt.png")
png(file = outfile, width=2800, height=2500, res=500)


cusHeat <- Heatmap(gg2@matrix, na_col = "white", col = gg2@matrix_color_mapping, bottom_annotation = gg2@bottom_annotation, left_annotation = gg2@left_annotation, 
                  show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = "left",row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                  column_title = "Differenital interaction strength", column_title_gp = gpar(fontsize = 10), column_names_rot = 60,
                  row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, #at = colorbar.break,
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
                 )

cusHeat
dev.off()


### Fig 8d/8e + Fig supp 8b-d: circos plots & violin plots
pathwayz <- c("FASLG","FN1","SEMA7")
pathwayz <- c("HGF","IGF","MPZ")
pathways <- pathwayz
subName <- ""

cellchat <- cellchat.oa
cellchat <- cellchat.norm

colz <- gg_color_hue(length(levels(cellchat@idents)))
names(colz) <- levels(cellchat@idents)

lapply(pathwayz, function(pathway){
    
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
    colz2 <- colz
    colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("./output/", outName, "/", subName, "/", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate(cellchat, layout = "chord", signaling = pathway, color.use = colz2, remove.isolate = F, big.gap = 5) #, group = groupzNames2
    dev.off()

    #get the active features in the pathway
    genez <- lapply(pathways, function(x){extractEnrichedLR(cellchat, signaling = x, geneLR.return = TRUE, enriched.only = T)[["geneLR"]]})
    names(genez) <- pathways

    Idents(seu.obj) <- "celltype.l2"
    seu.obj.sub <- subset(seu.obj, idents = cellTypeOFinterest)
    
    #plot expression using Seurat function
    pi <- VlnPlot(
        object = seu.obj.sub,
        pt.size = 0,
        same.y.lims = T,
        flip = T,
        group.by = "celltype.l2",
        fill.by = "ident",
        split.by = "cellSource",
#         cols = colz2,
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "/", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
})