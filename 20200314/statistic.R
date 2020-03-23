library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")
parser$add_argument("--seurat",
		    type="character")

args <- parser$parse_args()
if(!dir.exists(args$seurat)){
        dir.create(args$seurat,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
Idents(seurat)<-seurat$status
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1)






