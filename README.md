# equine_synovialFluid_scRNA

[![DOI](https://zenodo.org/badge/635774225.svg)](https://zenodo.org/doi/10.5281/zenodo.10680754)


This GitHub repository contains all the analysis code used in, "Characterization of the single cell landscape in normal and osteoarthritic equine joints."

The manuscript has been conditionally accepted  and will be published in Annals of Translational Medicine shortly.

## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the corresponding author(s), and/or contact Dylan Ammons at dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: input](/scripts/input) instructions for obtaining data associated with this study
- [:file\_folder: scripts](/scripts) contains the analysis code and source file used to complete the data analysis
- [:file\_folder: output](/scripts/output) contains the expected directory structure


## Supplemental data and potential uses:
0. [Browse the data](#browse-the-complete-annotated-dataset)
0. [Cell type annotations](#cell-type-annotations-with-defining-markers)
0. [Reference Mapping](#using-the-data-to-complete-reference-mapping)
0. [Module scoring](#module-scoring)

### Browse the complete annotated dataset

The proccessed dataset is avaliable for browsing via the UCSC Cell Browser portal.
Using the portal you can explore feature expression throughout the dataset as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Link to the dataset: https://equine-joint-oa.cells.ucsc.edu

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

### 2. Cell type annotations with defining markers

Cell markers lists were curated using the full dataset of 3 healthy and 4 chronic inflammatory enteropathy (CIE) dogs. The top 50 defining features (identified using `FindAllMarkers` for each cell type were considered, with the top 24 features evaluated for specificity using violin plots and preference given to unique features only found in the top 50 of one cell type.

<details open><summary>Cell type gene signatures</summary>
<p>

|Cell type     |                |Major markers                                                                                        |
|--------------|----------------|-----------------------------------------------------------------------------------------------|
|Macrophage    |                |APOE, PLTP, C1QC, CD68, ENSECAG00000024790, LGMN                                               |
|              |CCL2_Macrophage |CCL2, CCL8, LYVE1, P2RY13, HBEGF, PCP4L1                                                       |
|              |CD5L_Macrophage |MARCO, C4BPA, C1QC, CD55, MSR1, C1QB, ENSECAG00000022247                                       |
|              |GPNMB_Macrophage|APOE, CTSD, GPNMB, LGMN, ENSECAG00000024790, PLD3, ACP5                                        |
|              |TPPP3_Macrophage|ANXA2, S100A10, S100A4, ANXA1, LGALS3, TPPP3, S100A11, ANXA5                                   |
|              |Mo-Mac          |LYZ, PLAC8B, ENSECAG00000000436, SLPI, SERPINB10, AQP9                                         |
|              |Monocyte        |S100A12, ENSECAG00000010117, ECATH-3, VCAN, ENSECAG00000010598, ENSECAG00000010615             |
|Dendritic cell|                |CPVL, CST3, DQA, ENSECAG00000039383, DRA, DRB, IRF8, NAPSA, FCER1A, FLT3                       |
|              |cDC1            |ENSECAG00000039383, ENSECAG00000039998, CPVL, IRF8, ENSECAG00000029928, GPIHBP1, DNASE1L3, FLT3|
|              |cDC2            |FCER1A, CD1A7, CD1E2, DRB, DRA, DQB, CD1C, CD1A1                                               |
|              |cycling_DC      |BRCA1, CDCA7, UHRF1, CLSPN, RRM2, BRCA2, HELLS, CDT1                                           |
|              |migDC           |RORC, CCR7, NGFR, CHN1, SLCO5A1, SLC38A1, MARCKSL1, IL4I1                                      |
|              |pDC             |GJB2, MS4A1, GJA1, FCRLA, OSBPL10, TCF4, IL18R1                                                |
|Neutrophil    |                |ENSECAG00000030548, SNX10, SOD2, ILT11B, IFIT2, CCRL2, S100P, ENSECAG00000036967, SELP         |
|CD4 T cells   |                |LTB, ENSECAG00000036014, ENSECAG00000000910, CD5, CD3E, CD2                                    |
|              |SELL_CD4        |SELL, KLF2, CCR7, LEF1, ENSECAG00000010622, ENSECAG00000040878, EEF1A1                         |
|              |Activated_CD4   |ENSECAG00000000910, ENSECAG00000036014, CALY, CD40LG, ICOS                                     |
|              |Treg            |ARID5B, FGL2, RALA, FOXO1, TNFRSF1B, CGA                                                       |
|CD8 T cells   |                |CCL5, GZMA, CD7, CTSW, NKG7, CD3E, CD2, CD3G, FASLG                                            |
|              |Effector_CD8    |GZMA, CTSW, CCL5, NKG7, DRB, KLRK1                                                             |
|              |DAPL1_CD8       |DAPL1, GZMK, GZMM, IDO1, IRF7                                                                  |
|              |CX3CR1_T        |CX3CR1, ENSECAG00000031322, ENSECAG00000022193, SCML4, ZEB2                                    |
|              |IFN-T           |MX1, IRF7, ISG20, ENSECAG00000033029, ISG15, XAF1, SAMD9L                                      |
|gd T cells    |                |KLRB1, GNLY, ENSECAG00000038560, TRDC, BLK, KLRF1, RORA, NCR3                                  |
|              |GZMA_gd_T       |ENSECAG00000038560, BLK, TRDC, ENSECAG00000039473, KLRB1, GZMA, ENSECAG00000014832, CTSW       |
|              |IL23R_gd_T1     |KLRF1, SCART1, TRDC, KLRB1, TNFRSF25, IL23R, TNFSF13B, RORA, IKZF2                             |
|              |IL23R_gd_T2     |ENSECAG00000038560, BLK, SCART1, RHEX, KLRF1, TRDC, TNFRSF25, IL23R, RORA                      |
|              |KLRD1_T         |KLRB1, HOPX, CRYBG2, KLRD1, CCL5, KLRK1, TRDC, CTSW                                            |
|              |GNLY_T          |GNLY, KLRB1, KLRF1, CTSW, NKG7, TRDC, LAG3                                                     |
|B cells       |                |ENSECAG00000039599, MS4A1, CD79A, TCF4, CD19, PRAG1                                            |
|Cycling cells |                |PAFAH2, TOP2A, ENSECAG00000036105, H2AZ1, TUBA1A, NUSAP1, CENPF, H1-3                          |


</p>
</details>

### Using the data to complete reference mapping
Reference mapping is useful tool to facilitate the identification of cell types in single cell datasets. The approach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the reference datasets generated in this study. The default approach describes how to use the healthy only dataset, but it will also work with the combined dataset if you load that file in as the reference.

Before running the reference mapping code, a Seurat object need to be preprocessed and stored as an object named `seu.obj`.
```r
#set the path to the location in which the reference file is saved
reference <- readRDS(file = "eq_synovial_fluid_annotated.rds")

#prepare the reference
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#find conserved anchors with query and reference
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#select meta.data slot to use for label transfer -- change refdata value to use alternate labels (i.e., refdata = reference$celltype.l1)
predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$celltype.l3,
    dims = 1:50
)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

#generate and save the image
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "predicted.id",
    pt.size = 0.25,
    label = T,
    label.box = T,
    shuffle = F
)
ggsave("./output/referenceMap.png", width = 7, height = 7)
```

### 4. Module scoring
Module scoring is a supplemental approach that can be applied to single cell datasets with the goal of providing further insights into cell identities. The approach described below uses the Seurat function `AddModuleScore` and the gene lists presented in Table 2 of our associated manuscript. 

The concept of the AddModuleScore() function is similar to GSEA, but also distinct in many ways. Read the [Seurat documentation](https://satijalab.org/seurat/reference/addmodulescore) and/or check out [this webpage](https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/) for more details.

```r
#load in the reference file
# can download with this command: wget https://raw.githubusercontent.com/dyammons/equine_synovialFluid_scRNA/main/scripts/input/genesig_long.csv
ref.df <- read.csv("genesig_long.csv", header = T)

#organize the data
modulez <- split(ref.df$gene, ref.df$celltype.l2)

#complete module scoring
seu.obj <- AddModuleScore(seu.obj, features = modulez, name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results
features <- names(modulez)
ecScores <- DotPlot(
    seu.obj,
    assay = "RNA",
    features = features
)

ggsave("./output/dots_celltypes.png", width = 10, height = 6)
```

