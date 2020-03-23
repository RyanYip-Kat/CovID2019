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

Idents(seurat)<-seurat$UMAP_clusters
seurat<-RenameIdents(seurat,"9"="CDC1","12"="pDC","13"="CDC2","16"="pDC",
                     "2"="CD14","3"="CD14","4"="CD14","6"="CD14","7"="CD14","14"="CD14",
                     "5"="CD16","15"="CD16","17"="CD16",
                     "1"="Med","8"="Med","10"="Med","11"="Med")

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
	color = colorRampPalette(c("navy","#40E0D0","#B0171F"))(50)
        p<-pheatmap(cts,color=color,show_colnames =F,show_rownames = T,
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
