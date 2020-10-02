plot_figure=function(gene,dataset){
	setwd("/Users/TK/LRZ Sync+Share/Share_Bioinfo_Korn (Gildas Lepennetier)/2020.01.28_ScRNAseq_i-stream_vs_m-stream_3prime-5prime-combined")
	require(Seurat)
	require(ggplot2)
	#"3-prime" # OR	"5-prime"	"combined5+3"
	if(dataset=="3-prime"){GEX=readRDS("GEX_3.Rds");NAMETAG="_3prime"}
	if(dataset=="5-prime"){GEX=readRDS("GEX_5.Rds");NAMETAG="_5prime"}
	if(dataset=="combined5+3"){GEX=readRDS("GEX.singlet.integrated.rds");NAMETAG="_Combined5+3"}
	FEATURE=gene
	PLOT=FeaturePlot(GEX,features=FEATURE,pt.size=3,order=TRUE) + 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 18), #this is for the legend
			  axis.text = element_text(size = 26),
			  axis.ticks = element_line(size=1))
	OUTNAME=paste0("plots_Expression/",FEATURE,NAMETAG,".pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off();print(paste("saved:",OUTNAME))
	Idents(GEX) <- "stream"
	PLOT=VlnPlot(GEX,features=FEATURE,log=TRUE,pt.size = 0) + #slot = "counts", 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 26),
			  axis.text = element_text(size = 26),
			  axis.title.x = element_blank(),
			  axis.ticks = element_line(size=1),
			  legend.position = 'none') #+ scale_y_continuous(limits = c(0, NA))
	OUTNAME=paste0("plots_Expression/",FEATURE,NAMETAG,"_Vplot_stream.pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off();print(paste("saved:",OUTNAME))
	Idents(GEX) <- "group"
	PLOT=VlnPlot(GEX,features=FEATURE,log=TRUE,pt.size = 0) + #slot = "counts", 
		theme(axis.line = element_line(size=1),
			  text = element_text(size = 26),
			  axis.text = element_text(size = 26),
			  axis.title.x = element_blank(),
			  axis.ticks = element_line(size=1),
			  legend.position = 'none') #+ scale_y_continuous(limits = c(0, NA))
	OUTNAME=paste0("plots_Expression/",FEATURE,NAMETAG,"_Vplot_group.pdf");pdf(OUTNAME,width=10,height=10);print(PLOT);dev.off();print(paste("saved:",OUTNAME))
	
}
