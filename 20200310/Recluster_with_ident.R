library(Seurat)
library(monocle3)
library(dplyr)
library(argparse)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("-k",
                      type="integer",
                      default="20",
                      help="")

parser$add_argument("--ident",
                    type="integer",
                    nargs="+",
                    default=4,
                    help="idents to be recluster")


parser$add_argument("--outdir",
                     type="character",
                     default="")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/PCV_CAT_8/monocle3/model/20200309_seurat.rds")
cds<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/PCV_CAT_8/monocle3/model/monocle.rds")
Idents(seurat)=seurat$new_umap_clusters
new_cluster_list=list()
for(ident in args$ident){
	cells=WhichCells(seurat,idents=ident)
	cds_subset=cds[,cells]
	cds_subset=cluster_cells(cds_subset,
				 reduction_method="UMAP",
				 k=as.integer(args$k),
				 cluster_method="leiden",
                                 partition_qval=0.05)
	new_clusters=clusters(cds_subset,reduction_method="UMAP")
	new_clusters=as.data.frame(new_clusters)
	colnames(new_clusters)<-"new_cluster"
	x<-paste(ident,new_clusters$new_cluster,sep=".")
	new_clusters$new_cluster<-x

	new_cluster_list[[ident]]=new_clusters
}

cells=WhichCells(seurat,idents=args$ident,invert=TRUE)
meta=seurat@meta.data[cells,"new_umap_clusters",drop=FALSE]
colnames(meta)="new_cluster"
new.meta=do.call(rbind,new_cluster_list)

meta.data=rbind(meta,new.meta)
seurat<-AddMetaData(seurat,meta.data,col.name="new_cluster")
saveRDS(seurat,file.path(args$outdir,"20200310_seurat.rds"))
