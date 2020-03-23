suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

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
Idents(seurat)<-seurat$new_umap_clusters
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
                          min.pct = 0.2,
                          pseudocount.use = 1 )



if(!dir.exists(args$outdir)){
       dir.create(args$outdir,recursive=TRUE)
}
topn=markers%>%group_by(cluster)%>%top_n(70,wt=avg_logFC)
saveRDS(markers,file.path(args$outdir,"new_markers.rds"))
write.table(topn,file.path(args$outdir,"new_top70_markers.csv"),quote=F,row.names=F,sep=",")

