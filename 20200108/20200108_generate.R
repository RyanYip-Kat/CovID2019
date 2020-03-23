suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))


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
genes<-read.table("sample/su_genes.txt",stringsAsFactors=F)$V1
seurat<-readRDS(args$seurat)
seurat<-subset(seurat,CD3D==0 & CD3E==0 & CD3G==0 & CD19==0 & CD14==0 & CD7>0)

seurat<-FindNeighbors(seurat,features=genes)
seurat<-FindClusters(seurat,resolution=0.05)
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )



model.dir<-args$outdir
if(!dir.exists(model.dir)){
       dir.create(model.dir,recursive=TRUE)
}

pheatmap.idents<-function(seurat,features,slot="scale.data"){
        cts <- GetAssayData(seurat, slot = slot)
        if(slot=="counts"){
                cts <- log10(cts + 1)
        }
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(Idents(seurat))
        cts <- as.matrix(cts[features, names(new_cluster)])
        ac=data.frame(cluster=new_cluster)
        p<-pheatmap(cts,show_colnames =F,show_rownames = T,
                    cluster_rows = F,
                    cluster_cols = F,
                    fontsize=8,
                    annotation_col=ac)
        return(p)
}

topn<-markers%>%group_by(cluster)%>%top_n(15,wt=avg_logFC)
features<-as.character(topn$gene)
jpeg(file.path(args$outdir,"scale_doheatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=features,label=T)
dev.off()

jpeg(file.path(args$outdir,"counts_doheatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=features,label=T,slot="data")
dev.off()

jpeg(file.path(args$outdir,"scale_pheatmap.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,features)
dev.off()

jpeg(file.path(args$outdir,"counts_pheatmap.jpeg"),width=1024,height=1024)
pheatmap.idents(seurat,features,"data")
dev.off()
