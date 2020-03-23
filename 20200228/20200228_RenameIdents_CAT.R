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


seurat=readRDS("all_gene/CAT_8/monocle3/model/seurat.rds")
Idents(seurat)<-seurat$UMAP_clusters

seurat<-RenameIdents(seurat,"1"="Monocyte",
                     "2"="Center memory CD4",
                     "3"="Naive CD4 CD8",
                     "4"="CTL CD8/CTL CD4",
                     "5"="NK",
                     "6"="NK",
                     "7"="Transition CD4",
                     "8"="Effector memory CD8",
                     "9"="Naive mature B",
                     "10"="Monocyte",
                     "11"="Activated memory B",
                     "12"="DC4",
                     "13"="Monocyte",
                     "14"="Monocyte",
                     "15"="CTL CD4",
                     "16"="CD1C+ cDC2",
                     "17"="Effector memory CD4",
                     "18"="Erythroblasts",
                     "19"="ILC",
                     "20"="Megakaryocyte",
                     "21"="STMN1 + T",
                     "22"="Pulp cells",
                     "23"="CSF3R-hi CD14-hi Monocyte",
                     "24"="pDC",
                     "25"="NK",
                     "26"="Basophil",
                     "27"="CTL CD8",
                     "28"="Naive CD8",
                     "29"="STMN1 + T",
                     "30"="Immature B Cell",
                     "31"="Naive CD4",
                     "32"="cDC1",
                     "33"="NK",
                     "34"="Pulp cells"
                     )

saveRDS(seurat,"all_gene/CAT_8/monocle3/model/seurat.rds")
pdf(file.path(args$outdir,"CellType.pdf"))
DimPlot(seurat,reduction="umap",label=T)+theme(legend.text=element_text(size=7))
dev.off()


