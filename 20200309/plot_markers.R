library(Seurat)
library(ggplot2)
library(ggpubr)
library(argparse)


print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")



parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)

genes<-c("CD4","CD8A", "CCR7", "IL7R", "FOXP3","CXCR3",
	 "CCR6", "CCR4","GZMB","GZMK","RORC", "IL2RA", "CTLA4",
	 "PDCD1","CD274","ENTPD1","CD38","NT5E","CD24", "CD27","PF4")

for(gene in genes){
	p<-FeaturePlot(seurat,features=gene,reduction="tsne")
        jpeg(file.path(args$outdir,paste0(gene,"_tsne.jpeg")),width=1024,height=1024)
	print(p)
	dev.off()

	p<-FeaturePlot(seurat,features=gene,reduction="umap")
	jpeg(file.path(args$outdir,paste0(gene,"_umap.jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}


