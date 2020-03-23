library(Seurat)
library(argparse)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="seurat path")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}
seurat<-readRDS(args$seurat)

jpeg(file.path(args$outdir,"umap_clusters.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,
           reduction="umap",
           group.by="UMAP_clusters",
           label=T,label.size=7.5)
dev.off()

jpeg(file.path(args$outdir,"tsne_clusters.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,
             reduction="tsne",
             group.by="UMAP_clusters",
             label=T,label.size=7.5)
dev.off()

jpeg(file.path(args$outdir,"tsne_status.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,
	reduction="tsne",
	group.by="status")
dev.off()

jpeg(file.path(args$outdir,"umap_status.jpeg"),width = 1024,height = 1024)
DimPlot(seurat,
	reduction="umap",
	group.by="status")
dev.off()
