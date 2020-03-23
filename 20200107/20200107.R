library(stringr)
library(dplyr)
library(Seurat)
library(pheatmap)
library(argparse)
library(ggplot2)

print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the dataset  to be used")


parser$add_argument("--marker",
                    type="character",
                    default="",
                    help="the marker  to be used")


parser$add_argument("--outdir",
                      type="character",
                      default="",
                      help="save path")

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

pheatmap.idents<-function(seurat,features,slot="scale.data"){
	cts <- GetAssayData(seurat, slot = slot)
	if(slot=="counts"){
		cts <- log10(cts + 1)
	}
	features<-features[features%in%rownames(cts)]
	new_cluster <- sort(Idents(seurat))
	cts <- as.matrix(cts[features, names(new_cluster)])
	ac=data.frame(cluster=new_cluster)
	p<-pheatmap(cts,show_colnames =F,show_rownames = T,
		    cluster_rows = F,
		    cluster_cols = F,
		    fontsize=8,
		    annotation_col=ac)
	return(p)
}




seurat<-readRDS(args$seurat)
nk_genes<-read.table("sample/su_genes.txt",stringsAsFactors=F)$V1
markers<-readRDS(args$marker)
topn<-markers%>%group_by(cluster)%>%top_n(9,wt=avg_logFC)
markers<-as.character(topn$gene)
features<-setdiff(union(markers,nk_genes),c("RPS10","RPS17","RPS29","RPS20","RPL27A","RPL13A"))


jpeg(file.path(args$outdir,"scale_doheatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=features,label=T)
dev.off()

jpeg(file.path(args$outdir,"counts_doheatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=features,label=T,slot="data")
dev.off()





