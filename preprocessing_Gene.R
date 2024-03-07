
# Author:Songtao Wei
# Name:preprocessing_Gene

library(dplyr)
library(Seurat)
library(patchwork)

dataFileName='Gene_data_afterraw.txt';
GE_QC.data<- read.delim(dataFileName,header = TRUE,sep = "\t",row.names ="geneNames", check.names=FALSE)
GE_QC <- CreateSeuratObject(counts = GE_QC.data, min.cells = 5, min.features = 50);
GE_QC
VlnPlot(GE_QC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot1 <- FeatureScatter(GE_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GE_QC <- subset(GE_QC, subset = nFeature_RNA > 50 & nFeature_RNA < 4000)
GE_QC <- NormalizeData(GE_QC, normalization.method = "LogNormalize", scale.factor = 10000)



GE_QC <- FindVariableFeatures(GE_QC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GE_QC)
GE_QC <- ScaleData(GE_QC, features = all.genes)
Out_data <- as.data.frame(as.matrix(GE_QC@assays$RNA@scale.data))
write.csv(Out_data, file = "Gene_after_scaling.csv")
Out_Gene <- as.matrix(as.matrix(GE_QC@assays$RNA@var.features))
write.csv(Out_Gene, file = "Gene_List_selected.csv")
#all.genes <- rownames(GE_sel)


#GE_sel <- ScaleData(GE_sel, features = GE_sel)
#tmp <- GE_sel@assays[["RNA"]]@data
#tmp <- as.matrix(tmp)

