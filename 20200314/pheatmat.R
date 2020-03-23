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

seurat<-readRDS(args$seurat)
markers<-readRDS(args$marker)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
markers<-markers[!rownames(markers)%in%drop_genes,]
topn<-markers%>%group_by(cluster)%>%top_n(20,wt=avg_logFC)
genes<-c(topn$gene,"CD8A","CD14","CD68","CD19")


pheatmap.idents<-function(seurat,features,slot="scale.data"){
        cts <- GetAssayData(seurat, slot = slot)
        if(slot=="counts"){
                cts <- log10(cts + 1)
        }
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(seurat$status)
        cts <- as.matrix(cts[features, names(new_cluster)])
        ac=data.frame(cluster=new_cluster)
	#color<-c("#40E0D0","#B0171F")
        p<-pheatmap(cts,show_colnames =F,show_rownames = T,
                    cluster_rows = F,
                    cluster_cols = F,
                    fontsize=8,
                    annotation_col=ac)
        return(p)
}


jpeg(file.path(args$outdir,"heatmap1.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes)
dev.off()

jpeg(file.path(args$outdir,"heatmap2.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes,"data")
dev.off()
