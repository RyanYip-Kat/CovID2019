library(Seurat)
library(argparse)
print("...conplot parameters...")
parser <- ArgumentParser(description='Process some tasks')
#parser$add_argument("--markers",
#                    type="character",
#                    nargs="+",
#                    default=NULL,
#                    help="markers to be featureplot")

parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="which sample txt file to be used")

parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


seurat<-readRDS(args$seurat)
x<-subset(seurat,CD3E>0 & CD4>0)

genes.1<-c("JAK1","STAT3")
genes.2<-c("IL2RG","IL6ST","CXCR4")

jpeg(file.path(args$outdir,"JAK1_STAT3.jpeg"),width=1024,height=1024)
VlnPlot(x,features=genes.1,group.by="sample",ncol=1)
dev.off()

jpeg(file.path(args$outdir,"IL2RG_IL6ST_CXCR4.jpeg"),width=1024,height=1024)
VlnPlot(x,features=genes.2,group.by="sample",ncol=1)
dev.off()

