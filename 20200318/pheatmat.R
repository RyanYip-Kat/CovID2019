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
seurat$status<-factor(seurat$status,levels=c("HC","ER","TR"))
markers<-readRDS(args$marker)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
markers<-markers[!rownames(markers)%in%drop_genes,]
topn<-markers%>%group_by(cluster)%>%top_n(30,wt=avg_logFC)
genes<-topn$gene


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
	color = colorRampPalette(c("navy","#40E0D0","#B0171F"))(50)
        p<-pheatmap(cts,color=color,show_colnames =F,show_rownames = T,
                    cluster_rows = F,
                    cluster_cols = F,
                    fontsize=10,
                    annotation_col=ac)
        return(p)
}


pdf(file.path(args$outdir,"heatmap.pdf"),width=8,height=12)
pheatmap.idents(seurat,genes,"data")
dev.off()

