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
cells<-colnames(seurat)
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
status<-ifelse(ident%in%c(1,2,3,4),"Normal",ifelse(ident%in%c(12,9,6,8),"LT","ST"))
seurat$status<-status

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
Idents(seurat)<-seurat$status
seurat<-subset(seurat,idents=c("Normal","ST"))
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
			  features = keep_genes,
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )



if(!dir.exists(args$outdir)){
       dir.create(args$outdir,recursive=TRUE)
}
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
pheatmap.idents(seurat,genes)
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap2.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,genes,"data")
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap3.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes)
dev.off()

jpeg(file.path(args$outdir,"NST_heatmap4.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,slot="data")
dev.off()
