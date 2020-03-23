suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
library(pheatmap)
library(ggplot2)
print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")


args<-parser$parse_args()
if(!dir.exists(args$outdir)){
         dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat<-readRDS("../20200315/output/MC/model/seurat.rds")

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$UMAP_clusters
seurat<-subset(seurat,idents=c(9,12,13,16))
Idents(seurat)<-seurat$status
x<-subset(seurat,idents=c("HC","ER"))
#levels(x)<-c("ER","HC")
x$status<-factor(x$status,levels=c("ER","HC"))
markers <- FindAllMarkers(x, only.pos = FALSE,
                            test.use = "t",
                            features = keep_genes,
                            min.pct = 0.2,
                            logfc.threshold = 0.1,
                            pseudocount.use = 1 )


print(head(markers))
write.table(markers,file.path(args$outdir,"ER_HC_markers.csv"),quote=F,row.names=F,sep=",")

topn=markers%>%group_by(cluster)%>%top_n(25,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"HC_ER_heatmap1.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,size=7)+theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()

jpeg(file.path(args$outdir,"HC_ER_heatmap2.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,slot="data",size=7)+theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()

x<-subset(seurat,idents=c("ER","TR"))
#levels(x)<-c("ER","TR")
x$status<-factor(x$status,levels=c("ER","TR"))
markers <- FindAllMarkers(x, only.pos = FALSE,
                              test.use = "wilcox",
                              features = keep_genes,
                              min.pct = 0.2,
                              logfc.threshold = 0.1,
                              pseudocount.use = 1 )



write.table(markers,file.path(args$outdir,"ER_TR_markers.csv"),quote=F,row.names=F,sep=",")
topn=markers%>%group_by(cluster)%>%top_n(25,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"ER_TR_heatmap1.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,size=7)+
        theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()


jpeg(file.path(args$outdir,"ER_TR_heatmap2.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,slot="data",size=7)+theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()
