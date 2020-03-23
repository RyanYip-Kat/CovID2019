suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
library(pheatmap)

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
print("### Loading dataset")
seurat<-readRDS(args$seurat)
genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
#Idents(seurat)<-seurat$status
#seurat<-subset(seurat,idents=c("Normal","ST"))

Idents(seurat)<-seurat$UMAP_clusters
seurat<-RenameIdents(seurat,"4"="Naive CD4","5"="CD4 Tem","8"="CD4 Tem","9"="Naive CD4",
                     "10"="CD4 Tcm","12"="Treg","14"="CD4 Tcm","20"="Naive CD4","21"="Naive CD4",
                     "2"="CD8 CTL","6"="CD8 CTL","7"="CD8 Tm","11"="CD8 Tm","13"="Naive CD8",
                     "15"="CD8 CTL","16"="CD8 CTL","17"="CD8 CTL",
                     "1"="NK1","3"="NK1","19"="NK2",
                     "18"="Pro-T")
Idents(seurat)<-seurat$UMAP_clusters
seurat<-subset(seurat,idents=c(4,5,8,9,10,12,14,20,21))
Idents(seurat)<-seurat$status
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
			  features = keep_genes,
                          min.pct = 0.2,
                          logfc.threshold = 0.1,
                          pseudocount.use = 1 )



if(!dir.exists(args$outdir)){
       dir.create(args$outdir,recursive=TRUE)
}
topn=markers%>%group_by(cluster)%>%top_n(10,wt=avg_logFC)
saveRDS(markers,file.path(args$outdir,"celltype_markers.rds"))
write.table(markers,file.path(args$outdir,"celltype_markers.csv"),quote=F,row.names=F,sep=",")

pheatmap.idents<-function(seurat,features,slot="scale.data"){
        cts <- GetAssayData(seurat, slot = slot)
        if(slot=="counts"){
                cts <- log10(cts + 1)
        }
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(Idents(seurat))
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

genes<-topn$gene
jpeg(file.path(args$outdir,"celltype_heatmap1.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes)
dev.off()

jpeg(file.path(args$outdir,"celltype_heatmap2.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes,"data")
dev.off()

jpeg(file.path(args$outdir,"celltype_heatmap3.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes)
dev.off()

jpeg(file.path(args$outdir,"celltype_heatmap4.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,slot="data")
dev.off()
