library(Seurat)
library(ggplot2)
library(ggpubr)
library(argparse)


print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--s1",
                    type="character",
                    default="")

parser$add_argument("--s2",
                      type="character",
                      default="")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--n1",
                    type="character",
                    default="",
                    help="upgrated name")

parser$add_argument("--n2",
                      type="character",
                      default="",
                      help="upgrated name")
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat1<-readRDS(args$s1)
seurat2=readRDS(args$s2)

genes<-c("JAK1","JAK2","JAK3","TYK2",
	 "STAT3","STAT4", "STAT1","STAT2",
	 "STAT5A", "STAT6", "IL2RG","CXCR4",
	 "CXCL13", "CXCL12", "CCR10", "IL6ST","MKI67")


for(gene in genes){
	p1<-FeaturePlot(seurat1,features=gene,reduction="tsne")+ggtitle(paste0(gene,"(",args$n1,")"))
	p2<-FeaturePlot(seurat2,features=gene,reduction="tsne")+ggtitle(paste0(gene,"(",args$n2,")"))
        jpeg(file.path(args$outdir,paste0(gene,"_tsne.jpeg")),width=1024,height=1024)
	print(ggarrange(p1, p2,ncol=1, common.legend = TRUE))
	dev.off()

	p1<-FeaturePlot(seurat1,features=gene,reduction="umap")+ggtitle(paste0(gene,"(",args$n1,")"))
	p2<-FeaturePlot(seurat2,features=gene,reduction="umap")+ggtitle(paste0(gene,"(",args$n2,")"))
	jpeg(file.path(args$outdir,paste0(gene,"_umap.jpeg")),width=1024,height=1024)
	print(ggarrange(p1, p2,ncol=1, common.legend = TRUE))
	dev.off()
}


