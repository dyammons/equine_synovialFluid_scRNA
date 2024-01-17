#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(CellChat)

##############################
### CELLCHAT PREPROCESSING ###
##############################

#load in the annotated dataset
seu.obj <- readRDS("./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds")
seu.obj.backup <- seu.obj

### Run cell chat on OA samples
seu.obj <- subset(seu.obj, 
                  subset = cellSource == "OA")

#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                     gene_input = "rownames", 
                                     gene_output = "rownames", 
                                     input_species = "ecaballus",
                                     output_species = "human",
                                     non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/Nov_30_2023_cellChatobj_oa_ctl2.rds")


### Run cell chat on normal samples
seu.obj <- seu.obj.backup

Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Normal")

#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "ecaballus",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human 

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# saveRDS(cellchat, "./output/cellchat/Nov_30_2023_cellChatobj_norm_ctl2.rds")
saveRDS(cellchat, "./output/cellchat/Jan_16_2023_cellChatobj_norm_ctl2.rds") #updated with min.cells == 5



##################################
### END CELLCHAT PREPROCESSING ###
##################################


###############################
### BEGIN CELLCHAT ANALYSIS ###
###############################

#set output specifications
outName <- "cellchat"
subName <- "sf"

#restore full Seurat object
seu.obj <- seu.obj.backup
# seu.obj <- readRDS("./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds") #alterantively reload
# seu.obj.backup <- seu.obj

#load in processed cellchat data
# cellchat.norm <- readRDS("./output/cellchat/Nov_30_2023_cellChatobj_norm_ctl2.rds")
cellchat.norm <- readRDS("./output/cellchat/Jan_16_2023_cellChatobj_norm_ctl2.rds")
cellchat.oa <- readRDS("./output/cellchat/Nov_30_2023_cellChatobj_oa_ctl2.rds")
tc.df <- read.csv("sf_idents_tcell_12-05-2023.csv")
mye.df <- read.csv("sf_idents_myeloid_12-05-2023.csv")

#load in the cell type colors
colz.base <- c(tc.df$newCol, mye.df$newCol, gg_color_hue(8)[c(8,1)])
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, unique(seu.obj$celltype.l1)[c(9,8)])

#get overlapping pathways
print(cellchat.oa@netP$pathways[cellchat.oa@netP$pathways %in% cellchat.norm@netP$pathways])

#merge cellchat objects
object.list <- list(Normal = cellchat.norm, OA = cellchat.oa)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP")
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


### Fig supp 4a - bar chart with number of interactions
outfile <- paste0("./output/", outName, "/", subName, "_supp4a.png")
data.df <- compareInteractions(cellchat, show.legend = F)$data

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
ggsave(outfile, height = 1, width = 6)


### Fig 4a - interactivity scater plots by condition

#create plot using CellChat built in function
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
gg.df <- gg.df %>% group_by(data_type) %>% arrange(desc(strength)) %>% mutate(lab=ifelse(row_number() <= 5 | labels == "Monocyte", as.character(labels), NA)) %>% ungroup()

#make figure
pis <- lapply(c("Normal","OA"),function(z){
    gg.df.sub <- gg.df %>% filter(data_type == z)

    ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, colour = labels, label=lab)) + 
            ggtitle(z) +
            geom_point(color = gg.df.sub$colz.base) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black", box.padding = 1,  min.segment.length = 0.25, max.iter = 1000000, seed = 666) + 
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
ggsave(pi2, file = paste0("./output/", outName, "/", subName, "_fig4a.png"), width = 3.25, height = 6)


### Fig 4b - create differential interaction heatmap

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

outfile <- paste0("./output/", outName, "/", subName, "_fig4b.png")
png(file = outfile, width=2800, height=2500, res=500)

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


### Fig supp - plot myeloid-myeloid specific chaanges

### Run cell chat on oa samples
seu.obj <- seu.obj.backup

### Run cell chat on OA samples
seu.obj <- subset(seu.obj, 
                  subset = cellSource == "OA")

Idents(seu.obj) <- "celltype.l1"
seu.obj <- subset(seu.obj, 
                 subset = celltype.l1 == "Macrophage" | celltype.l1 == "Monocyte" | celltype.l1 == "Neutrophil" | celltype.l1 == "DC")
# seu.obj <- subset(seu.obj, subset = celltype.l2 != "migDC")
#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "ecaballus",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human 

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/Jan_01_2024_cellChatobj_oa_MYELOID.rds")


### Run cell chat on normal samples
seu.obj <- seu.obj.backup

### Run cell chat on OA samples
seu.obj <- subset(seu.obj, 
                  subset = cellSource == "Normal")

Idents(seu.obj) <- "celltype.l1"
seu.obj <- subset(seu.obj, 
                 subset = celltype.l1 == "Macrophage" | celltype.l1 == "Monocyte" | celltype.l1 == "Neutrophil" | celltype.l1 == "DC")
# seu.obj <- subset(seu.obj, subset = celltype.l2 != "migDC")
#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "ecaballus",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l2")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human 

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/Jan_01_2024_cellChatobj_norm_MYELOID.rds")

#set output specifications
outName <- "cellchat"
subName <- "sf"

#restore full Seurat object
seu.obj <- seu.obj.backup
# seu.obj <- readRDS("./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds") #alterantively reload
# seu.obj.backup <- seu.obj

#load in processed cellchat data
cellchat.norm <- readRDS("./output/cellchat/Jan_01_2024_cellChatobj_norm_MYELOID.rds")
cellchat.oa <- readRDS("./output/cellchat/Jan_01_2024_cellChatobj_oa_MYELOID.rds")
tc.df <- read.csv("sf_idents_tcell_12-05-2023.csv")
mye.df <- read.csv("sf_idents_myeloid_12-05-2023.csv")

#load in the cell type colors
colz.base <- c(tc.df$newCol, mye.df$newCol, gg_color_hue(8)[c(8,1)])
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, unique(seu.obj$celltype.l1)[c(9,8)])
colz.base <- colz.base[13:24]
#get overlapping pathways
print(cellchat.oa@netP$pathways[cellchat.oa@netP$pathways %in% cellchat.norm@netP$pathways])

#merge cellchat objects
object.list <- list(Normal = cellchat.norm, OA = cellchat.oa)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP")
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#plot using CellChat
i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)


#extract the data and customize the heatmap
colz.base <- colz.base[match(colnames(ht1@matrix), names(colz.base))]

ha <- HeatmapAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )

cusHeat1 <- Heatmap(ht1@matrix, na_col = "white", col = ht1@matrix_color_mapping, 
                   bottom_annotation = ha,
                   show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                   column_title = "Outgoing signal strength (Normal)", column_title_gp = gpar(fontsize = 10), column_names_rot = 60,
                   row_title = "Network", row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                   cluster_rows=F, cluster_columns=F,
                  
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                               border = NA,
                                               legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
                 )

colz.base <- c(tc.df$newCol, mye.df$newCol, gg_color_hue(8)[c(8,1)])
names(colz.base) <- c(tc.df$celltype.l2, mye.df$majorID, unique(seu.obj$celltype.l1)[c(9,8)])
colz.base <- colz.base[13:24]

colz.base <- colz.base[match(colnames(ht2@matrix), names(colz.base))]

ha <- HeatmapAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )

cusHeat2 <- Heatmap(ht2@matrix, na_col = "white", col = ht2@matrix_color_mapping, 
                   bottom_annotation = ha,
                   show_column_dend = FALSE, show_row_dend = FALSE, row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                   column_title = "Outgoing signal strength (OA)", column_title_gp = gpar(fontsize = 10), column_names_rot = 60,
                   row_title = "Network", row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                   cluster_rows=F, cluster_columns=F,
                  
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                               border = NA,
                                               legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
                 )

png(file = paste0("./output/", outName, "/", subName, "_", i, "_test.png"), width=3000, height=2000, res=500)
draw(cusHeat1 + cusHeat2, ht_gap = unit(0.5, "cm"))
dev.off()


### Fig supp 4b - unstacked bargraph of network information flow

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
                 # plot.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key.size = unit(1,"line"),
                  axis.title = element_text(size = 18),
                  axis.text = element_text(size = 16),
                  plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
            ) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0,0)) + 
            labs(y = "Interaction network",
                 x = "Relative information flow",
                 title = "Myeloid-myeloid interactions") + 
            guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + NoLegend()

ggsave(p, file = paste0("./output/", outName, "/", subName, "_fig4c.png"))








### Fig 4c - stacked bargraph of network information flow

#plot information flow using CellChat built-in function
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

#extract data and plot using custom code
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

ggsave(p, file = paste0("./output/", outName, "/", subName, "_fig4c.png"))


### Fig supp 4b - unstacked bargraph of network information flow

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

ggsave(p, file = paste0("./output/", outName, "/", subName, "_supp4c.png"))


### Fig 4d/4e + Fig supp 4c/d - circos plots & violin plots
#fix annotation error
Idents(seu.obj) <- "celltype.l1"
seu.obj <- RenameIdents(seu.obj, c("Monocyte" = "Macrophage"))
seu.obj$celltype.l1 <- Idents(seu.obj)

groupz.df <- table(seu.obj$celltype.l1, seu.obj$celltype.l2) %>% melt() %>% filter(value > 0)
groupzNames <- groupz.df$Var.1
names(groupzNames) <- groupz.df$Var.2

cellchat <- cellchat.oa
pathwayz <- c("FASLG","JAM","VISFATIN")
pathways <- pathwayz
subName <- "sf"

cellchat <- cellchat.norm
pathwayz <- c("HGF","IGF")
pathways <- pathwayz
subName <- "sf"

pis <- lapply(pathwayz, function(pathway){
    
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
    outfile <- paste("./output/", outName, "/", subName, "_", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate_mod(cellchat, layout = "chord", signaling = pathway, color.use = colz2, remove.isolate = F, big.gap = 5, group = groupzNames, title.name = NULL)
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
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "_", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
    return(pi)
    
}) 

#manually clean the figures - oa
pis[[1]] + ggtitle(paste0(pathways[1], " network")) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1"))
ggsave(paste("./output/", outName, "/", subName, "_", pathways[1] ,"_viln.png", sep = ""), height = 3, width = 6)

pis[[2]] + ggtitle(paste0(pathways[2], " network")) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + theme(plot.margin = unit(c(7,7,7,21), "pt"))
ggsave(paste("./output/", outName, "/", subName, "_", pathways[2] ,"_viln.png", sep = ""), height = 3, width = 6)

pis[[3]] + ggtitle(paste0(pathways[3], " network")) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + theme(plot.margin = unit(c(7,7,7,49), "pt"))
ggsave(paste("./output/", outName, "/", subName, "_", pathways[3] ,"_viln.png", sep = ""), height = 4, width = 6)



#manually clean the figures - norm
pis[[1]] + ggtitle(paste0(pathways[1], " network")) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + theme(plot.margin = unit(c(7,7,7,21), "pt"))
ggsave(paste("./output/", outName, "/", subName, "_", pathways[1] ,"_viln.png", sep = ""), height = 3.5, width = 6)

pis[[2]] + ggtitle(paste0(pathways[2], " network")) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) 
ggsave(paste("./output/", outName, "/", subName, "_", pathways[2] ,"_viln.png", sep = ""), height = 3, width = 6)


#############################
### END CELLCHAT ANALYSIS ###
#############################



### Fig supp 4a - bar chart with number of interactions
outfile <- paste0("./output/", outName, "/", subName, "_supp4a.png")
data.df <- compareInteractions(cellchat, show.legend = F)$data

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
ggsave(outfile, height = 1, width = 6)


### Fig 4a - interactivity scater plots by condition

#create plot using CellChat built in function
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

#make figure
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
ggsave(pi2, file = paste0("./output/", outName, "/", subName, "_fig4a.png"), width = 3.25, height = 6)


### Fig 4b - create differential interaction heatmap

#create heatmap using CellChat built-in function
gg2 <- netVisual_heatmap(cellchat, measure = "count",
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

outfile <- paste0("./output/", outName, "/", subName, "_fig4b.png")
png(file = outfile, width=2800, height=2500, res=500)

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


### Fig 4c - stacked bargraph of network information flow

#plot information flow using CellChat built-in function
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

#extract data and plot using custom code
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

ggsave(p, file = paste0("./output/", outName, "/", subName, "_fig4c.png"))


### Fig supp 4b - unstacked bargraph of network information flow

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

ggsave(p, file = paste0("./output/", outName, "/", subName, "_fig4c.png"))


### Fig 4d/4e + Fig supp 4c/d - circos plots & violin plots
#fix annotation error
Idents(seu.obj) <- "celltype.l1"
seu.obj <- RenameIdents(seu.obj, c("Monocyte" = "Macrophage"))
seu.obj$celltype.l1 <- Idents(seu.obj)

groupz.df <- table(seu.obj$celltype.l1, seu.obj$celltype.l2) %>% melt() %>% filter(value > 0)
groupzNames <- groupz.df$Var.1
names(groupzNames) <- groupz.df$Var.2

cellchat <- cellchat.oa
pathwayz <- c("FASLG","FN1","SEMA7")
pathways <- pathwayz
subName <- "sf"

cellchat <- cellchat.norm
pathwayz <- c("HGF","IGF","MPZ")
pathways <- pathwayz
subName <- "sf"

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
    colz2 <- colz.base
    colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("./output/", outName, "/", subName, "_", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate_mod(cellchat, layout = "chord", signaling = pathway, color.use = colz2, remove.isolate = F, big.gap = 5, group = groupzNames, title.name = NULL)
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
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "_", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
})



lapply(c(1:12), function(i){
    png(file = paste0("./output/", outName, "/", subName, "_", i, "_test.png"))
    p <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:12), comparison = c(1, 2), remove.isolate = FALSE)
    print(p)
    dev.off()
})

