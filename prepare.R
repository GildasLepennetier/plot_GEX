#from: GEX.singlet.integrated.rds
#do:   GEX_5.Rds, GEX_3.Rds
require(Seurat)
require(ggplot2)
setwd("C:/Users/gilou/LRZ Sync+Share/Share_Bioinfo_Korn/2020.01.28_ScRNAseq_i-stream_vs_m-stream_3prime-5prime-combined")
GEX <- readRDS("GEX.singlet.integrated.rds")
group <- read.table("group_sc_info.txt", header = TRUE,row.names=1)
GEX <- AddMetaData(object=GEX, metadata=group)
DefaultAssay(object = GEX) <- "RNA"
GEX@meta.data$stream=substr(GEX@meta.data$group,1,3)
GEX@meta.data$origin = substr(GEX@meta.data$exp.group,1,1)
GEX_5 = subset(GEX, origin == "5")
GEX_3 = subset(GEX, origin == "3")
saveRDS(GEX_5,"GEX_5.Rds")
saveRDS(GEX_3,"GEX_3.Rds")
