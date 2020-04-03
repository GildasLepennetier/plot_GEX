plot_me_that=function(dir,SYMBOL){
	setwd(dir)
	require(Seurat)
	require(ggplot2)
	GEX.singlet.integrated <- readRDS("GEX.singlet.integrated.rds")
	group <- read.table("group_sc_info.txt", header = TRUE,row.names=1)
	GEX.singlet.integrated <- AddMetaData(object=GEX.singlet.integrated, metadata=group)
	DefaultAssay(object = GEX.singlet.integrated) <- "RNA"
	GEX.singlet.integrated@meta.data$stream=substr(GEX.singlet.integrated@meta.data$group,1,3)
	FEATURE=SYMBOL
	PLOT=FeaturePlot(GEX.singlet.integrated,features=FEATURE,pt.size=3,sort.cell=TRUE) + 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 18), #this is for the legend
			  axis.text = element_text(size = 26),
			  axis.ticks = element_line(size=1))
	OUTNAME=paste0(FEATURE,".pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off()
	Idents(GEX.singlet.integrated) <- "stream"
	PLOT=VlnPlot(GEX.singlet.integrated,features=FEATURE,log=TRUE,pt.size = 0) + #slot = "counts", 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 26),
			  axis.text = element_text(size = 26),
			  axis.title.x = element_blank(),
			  axis.ticks = element_line(size=1),
			  legend.position = 'none') #+ scale_y_continuous(limits = c(0, NA))
	OUTNAME=paste0(FEATURE,"_Vplot_stream.pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off()
	Idents(GEX.singlet.integrated) <- "group"
	PLOT=VlnPlot(GEX.singlet.integrated,features=FEATURE,log=TRUE,pt.size = 0) + #slot = "counts", 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 26),
			  axis.text = element_text(size = 26),
			  axis.title.x = element_blank(),
			  axis.ticks = element_line(size=1),
			  legend.position = 'none') #+ scale_y_continuous(limits = c(0, NA))
	OUTNAME=paste0(FEATURE,"_Vplot_group.pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off()
}
