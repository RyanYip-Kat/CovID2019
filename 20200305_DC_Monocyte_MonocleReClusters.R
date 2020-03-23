library(Seurat)
library(monocle3)
library(argparse)

print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--cds",
                      type="character",
                      default="",
                      help="the dataset  to be used")

parser$add_argument("--k",
		    type="integer",
		    default=20)

parser$add_argument("--outdir",
		    type="character",
		    default="",
		    help="the path to save object")

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


seurat<-readRDS(args$seurat)
Idents(seurat)=seurat$UMAP_clusters
#idents<-unique(c(2, 3, 4, 7, 15, 17, 31, 3, 8, 27, 28)) #CAT
#idents<-unique(c(1, 2, 4, 6, 8, 4, 6, 13, 14, 15)) # UN
#idents<-unique(c(1, 3, 4, 6, 9, 10, 21, 22)) # pcv
#idents<-unique(c(1, 2, 3, 8,10, 15, 2, 7, 13, 23)) # VHK
idents<-unique(c(3,5,6,14,30,2,4,10,12,13,16,17,25,28,31,32,35)) # BD
cells=WhichCells(seurat,idents=idents)
cds<-readRDS(args$cds)

cds_subset<-cds[,cells]
cds_subset<-cluster_cells(cds_subset,reduction_method="UMAP",cluster_method="leiden",partition_qval=0.01,k=as.integer(args$k))
UMAP_clusters<-clusters(cds_subset,reduction_method="UMAP")
print(table(UMAP_clusters))
UMAP_clusters<-as.data.frame(UMAP_clusters)
colnames(UMAP_clusters)<-"new_umap_clusters"


seurat<-subset(seurat,idents=idents)
seurat=AddMetaData(seurat,metadata=UMAP_clusters)
saveRDS(seurat,file.path(args$outdir,"seurat.rds"))


pdf(file.path(args$outdir,"umap_clusters.pdf"))
p1<-DimPlot(seurat,reduction="umap",group.by="UMAP_clusters",label=TRUE,label.size=7.5)
p2<-DimPlot(seurat,reduction="umap",group.by="new_umap_clusters",label=TRUE,label.size=7.5)
print(CombinePlots(list(p1,p2),ncol=1))
dev.off()
