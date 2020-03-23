library(Seurat)
library(argparse)
library(ggplot2)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


seurat=readRDS("all_gene/BD_8/monocle3/model/seurat.rds")
Idents(seurat)<-seurat$UMAP_clusters

seurat<-RenameIdents(seurat,"1"="NK",
                     "2"="CTL CD8",
                     "3"="Naive CD4",
                     "4"="Naive CD8",
                     "5"="Central memory CD4",
                     "6"="Central memory CD4",
                     "7"="Classical monocyte",
                     "8"="Naive B",
                     "9"="Memory B",
                     "10"="CTL CD8",
                     "11"="Monocyte",
                     "12"="Effector memory CD8",
                     "13"="Central memory CD8",
                     "14"="TCM CD4",
                     "15"="TNF-hi CD14-hi Monocyte",
                     "16"="Unknow",
                     "17"="KLRB1-hi effector memory CD8",
                     "18"="Monocyte",
                     "19"="Unknow",
                     "20"="NK",
                     "21"="DC4",
                     "22"="NK",
                     "23"="Memory B",
                     "24"="Monocyte",
                     "25"="CTL CD8",
                     "26"="IL1-hi TNF-hi Monocyte",
                     "27"="Megakaryocyte",
                     "28"="Unknow",
                     "29"="Unknow",
                     "30"="Unknow",
                     "31"="KLRB1-hi effector memory CD8",
                     "32"="Naive CD8",
                     "33"="Pulp cells",
                     "34"="CD1C+ cDC2",
                     "35"="Naive CD8",
                     "36"="Naive B",
                     "37"="Unknow")

saveRDS(seurat,"all_gene/BD_8/monocle3/model/seurat.rds")
pdf(file.path(args$outdir,"CellType.pdf"))
DimPlot(seurat,reduction="umap",label=T)+theme(legend.text=element_text(size=7))
dev.off()
