

seu.obj.backUp <- seu.obj


### run cell chat on treated
library(CellChat)
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "OA")

cnts <- seu.obj@assays$RNA@data
meta <- seu.obj@meta.data
meta$chatID <- meta$finalClusters
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "chatID") #  usin gpredicted ct.l2

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
saveRDS(cellchat, "./output/cellchat/cellChatobj_oa_finalClusters.rds")
# saveRDS(cellchat, "./output/cellchat/cellChatobj_oa_clusID.rds")


#run cell chat on transfered naive
seu.obj <- seu.obj.backUp

Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Normal")

table(seu.obj$cellSource)

cnts <- seu.obj@assays$RNA@data
meta <- seu.obj@meta.data
meta$chatID <- meta$finalClusters
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "chatID") #  using predicted ct.l2

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
saveRDS(cellchat, "./output/cellchat/cellChatobj_norm_finalClusters.rds")
# saveRDS(cellchat, "./output/cellchat/cellChatobj_norm_clusID.rds")



#### do the comparision analysis  --- this analysis may be confounded by differences in filter btwn tx and naive samples
outName <- "cellchat_syn"
subName <- ""

cellchat.norm <- readRDS("./output/cellchat/cellChatobj_norm_finalClusters.rds")
cellchat.oa <- readRDS("./output/cellchat/cellChatobj_oa_finalClusters.rds")

# cellchat.norm <- readRDS("./output/cellchat/cellChatobj_norm_clusID.rds")
# cellchat.oa <- readRDS("./output/cellchat/cellChatobj_oa_clusID.rds")

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
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + coord_flip()
#gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
#gg1 + gg2
ggsave(outfile, height=2,width=6)

#pretty worthless
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_net.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


#make heatmap -- pretty cool
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_weight.png")
png(file = outfile, width=2500, height=2500, res=400)
# gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T)
gg2
dev.off()

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_oaVn_heatMAP_weight.png")
png(file = outfile, width=2000, height=2200, res=500)


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

###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
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

# colz.df <- as.data.frame(colz)
# colz.df$labels <- rownames(colz.df)

# gg.df <- gg.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

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

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 2, guides = 'collect') & scale_size_continuous(limits = c(min(gg.df$Count), max(gg.df$Count))) & xlim(0, 12.5) & ylim(0, 12.5)
# pi <- pi1 + pi2 + plot_layout(ncol = 2, widths = c(0.22,0.78)) 
ggsave(pi2, file = paste0("./output/", outName, "/", outName, "_interactionScater.png"), width = 6.25, height = 3)







### Fig 8d/8e + Fig supp 8b-d: circos plots & violin plots
pathwayz <- c("CXCL")
pathwayz <- c("VCAM", "ICAM", "NOTCH")
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

    Idents(seu.obj) <- "finalClusters"
    seu.obj.sub <- subset(seu.obj, idents = cellTypeOFinterest)
    
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
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "/", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
})







