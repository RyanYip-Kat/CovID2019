suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
library(pheatmap)
library(ggplot2)
print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")


args<-parser$parse_args()
if(!dir.exists(args$outdir)){
         dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat<-readRDS(args$seurat)

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$UMAP_clusters
#seurat<-RenameIdents(seurat,"4"="CD4","5"="CD4","8"="CD4","9"="CD4",
#                     "10"="CD4","12"="CD4","14"="CD4","20"="CD4","21"="CD4",
#                     "2"="CD8","6"="CD8","7"="CD8","11"="CD8","13"="CD8",
#                     "15"="CD8","16"="CD8","17"="CD8",
#                     "1"="NK","3"="NK","19"="NK",
#                     "18"="Proliferation")

seurat<-subset(seurat,idents=c(2,3,4,5,6))
Idents(seurat)<-seurat$status
Idents(seurat)<-seurat$status
x<-subset(seurat,idents=c("HC","ER"))
levels(x)<-c("ER","HC")
#seurat$status<-factor(seurat$status,levels=c("HC","ER"))
markers <- FindAllMarkers(x, only.pos = FALSE,
                            test.use = "t",
                            features = keep_genes,
                            min.pct = 0.2,
                            logfc.threshold = 0.1,
                            pseudocount.use = 1 )



#saveRDS(markers,file.path(args$outdir,"CD14_HER_markers.rds"))
#write.table(markers,file.path(args$outdir,"CD14_HER_markers.csv"),quote=F,row.names=F,sep=",")

topn=markers%>%group_by(cluster)%>%top_n(25,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"HC_ER_heatmap.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,size=7)+theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()

x<-subset(seurat,idents=c("ER","TR"))
levels(x)<-c("ER","TR")
#seurat$status<-factor(seurat$status,levels=c("HC","ER"))
markers <- FindAllMarkers(x, only.pos = FALSE,
                              test.use = "t",
                              features = keep_genes,
                              min.pct = 0.2,
                              logfc.threshold = 0.1,
                              pseudocount.use = 1 )




topn=markers%>%group_by(cluster)%>%top_n(25,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"ER_TR_heatmap.jpeg"),width=1080,height=1024)
DoHeatmap(x,features=genes,size=7)+
        theme(axis.text.y=element_text(size=14),legend.text=element_text(size=12))
dev.off()
