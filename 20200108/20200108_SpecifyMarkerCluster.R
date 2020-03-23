suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparse))

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
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat<-readRDS(args$seurat)
seurat<-subset(seurat,CD3E==0 & CD3D==0 & CD3G==0 & CD14==0 & CD19==0 & CD7>0)
cells<-Cells(seurat)

################### SIGLEC10
x<-subset(seurat,NCAM1>1.0,slot="scale.data")
y<-subset(seurat,NCAM1<1.0,slot="scale.data")

seurat<-SetIdent(seurat,cells=Cells(x),value="NCAM1_HIGH")

genes<-read.table("sample/su_genes.txt",stringsAsFactors=F)$V1
y<-FindNeighbors(y,features=genes)
y<-FindClusters(y,resolution=0.05)

Idents(y)<-paste("NCAM1_LOW",Idents(y),sep="_")
idents<-unique(Idents(y))
for(ident in idents){
	cc<-WhichCells(y,idents=ident)
	seurat<-SetIdent(seurat,cells=cc,value=ident)
}
levels(seurat)<-levels(seurat)[order(levels(seurat))]
jpeg(file.path(args$outdir,"NCAM1_scale_heatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,label=FALSE)
dev.off()

jpeg(file.path(args$outdir,"NCAM1_counts_heatmap.jpeg"),width=1024,height=1024)
DoHeatmap(seurat,features=genes,slot="data",label=FALSE)
dev.off()
