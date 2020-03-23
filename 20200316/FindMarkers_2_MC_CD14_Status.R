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

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$UMAP_clusters
seurat<-RenameIdents(seurat,"9"="DC","12"="DC","13"="DC","16"="DC",
                     "2"="CD14","3"="CD14","4"="CD14","6"="CD14","7"="CD14","14"="CD14",
                     "5"="CD16","15"="CD16","17"="CD16",
                     "1"="Med","8"="Med","10"="Med","11"="Med")

seurat<-subset(seurat,idents="CD14")
Idents(seurat)<-seurat$status
seurat<-subset(seurat,idents=c("HC","ER"))
seurat$status<-factor(seurat$status,levels=c("HC","ER"))
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                            test.use = "t",
                            features = keep_genes,
                            min.pct = 0.2,
                            logfc.threshold = 0.1,
                            pseudocount.use = 1 )



saveRDS(markers,file.path(args$outdir,"CD14_HER_markers.rds"))
write.table(markers,file.path(args$outdir,"CD14_HER_markers.csv"),quote=F,row.names=F,sep=",")

topn=markers%>%group_by(cluster)%>%top_n(25,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"CD14_HCER_heatmap1.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes)
dev.off()

jpeg(file.path(args$outdir,"CD14_HCER_heatmap2.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,slot="data")
dev.off()

