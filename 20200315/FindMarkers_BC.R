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
if(!dir.exists(args$outdir)){
         dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat<-readRDS(args$seurat)
cells<-colnames(seurat)
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
status<-ifelse(ident%in%c(1,2,3,4),"Normal",ifelse(ident%in%c(12,9,6,8),"LT","ST"))
seurat$status<-status

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$UMAP_clusters
#seurat<-subset(seurat,idents=c(1,5,7)) # naive
#seurat<-subset(seurat,idents=c(2,4)) # Memory
#seurat<-subset(seurat,idents=c(3)) # Pb
seurat<-subset(seurat,idents=c(6)) # Immure
Idents(seurat)<-seurat$status
seurat1<-subset(seurat,idents=c("Normal","ST"))
markers <- FindAllMarkers(seurat1, only.pos = FALSE,
                          test.use = "t",
			  features = keep_genes,
                          min.pct = 0.2,
                          logfc.threshold = 0.1,
                          pseudocount.use = 1 )



topn=markers%>%group_by(cluster)%>%top_n(30,wt=avg_logFC)
saveRDS(markers,file.path(args$outdir,"NST_markers.rds"))
write.table(markers,file.path(args$outdir,"NST_markers.csv"),quote=F,row.names=F,sep=",")

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

genes<-topn$gene
jpeg(file.path(args$outdir,"NST_heatmap1.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat1,genes)
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap2.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat1,genes,"data")
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap3.jpeg"),width=1024,height=1024)
DoHeatmap(seurat1,features=genes)
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap4.jpeg"),width=1024,height=1024)
DoHeatmap(seurat1,features=genes,slot="data")
dev.off()


seurat2<-subset(seurat,idents=c("LT","ST"))
markers <- FindAllMarkers(seurat2, only.pos = FALSE,
                            test.use = "t",
                            features = keep_genes,
                            min.pct = 0.2,
                            logfc.threshold = 0.1,
                            pseudocount.use = 1 )



topn=markers%>%group_by(cluster)%>%top_n(30,wt=avg_logFC)
saveRDS(markers,file.path(args$outdir,"LST_markers.rds"))
write.table(markers,file.path(args$outdir,"LST_markers.csv"),quote=F,row.names=F,sep=",")

genes<-topn$gene
jpeg(file.path(args$outdir,"LST_heatmap1.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat2,genes)
dev.off()

jpeg(file.path(args$outdir,"LST_heatmap2.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat2,genes,"data")
dev.off()

jpeg(file.path(args$outdir,"LST_heatmap3.jpeg"),width=1024,height=1024)
DoHeatmap(seurat2,features=genes)
dev.off()

jpeg(file.path(args$outdir,"LST_heatmap4.jpeg"),width=1024,height=1024)
DoHeatmap(seurat2,features=genes,slot="data")
dev.off()

