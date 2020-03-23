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

df=read.csv("ERS_top3.csv")
print(dim(df))
barcode<-as.character(df$Barcode)

Idents(seurat)=seurat$status
seurat<-subset(seurat,idents="ER")
cells<-Cells(seurat)
celltype<-ifelse(cells%in%barcode,"ERS TOP3","Others")
print(table(celltype))
seurat$celltype<-celltype

Idents(seurat)<-celltype
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                            test.use = "t",
                            features = keep_genes,
                            min.pct = 0.2,
                            logfc.threshold = 0.1,
                            pseudocount.use = 1 )



saveRDS(markers,file.path(args$outdir,"markers.rds"))
write.table(markers,file.path(args$outdir,"markers.csv"),quote=F,row.names=F,sep=",")

topn=markers%>%group_by(cluster)%>%top_n(50,wt=avg_logFC)
genes<-topn$gene
jpeg(file.path(args$outdir,"heatmap1.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes)
dev.off()

jpeg(file.path(args$outdir,"heatmap2.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,slot="data")
dev.off()


pheatmap.idents<-function(seurat,features,slot="scale.data"){
        cts <- GetAssayData(seurat, slot = slot)
        if(slot=="counts"){
                cts <- log10(cts + 1)
        }
        features<-features[features%in%rownames(cts)]
        #new_cluster <- sort(seurat$status)
        new_cluster <- sort(seurat$celltype)
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


jpeg(file.path(args$outdir,"heatmap3.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes)
dev.off()

jpeg(file.path(args$outdir,"heatmap4.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes,"data")
dev.off()

