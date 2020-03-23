suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
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
if(!dir.exists(args$outdir)){
	dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat<-readRDS(args$seurat)
cells<-colnames(seurat)
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
status<-ifelse(ident%in%c(1,2,3,4),"Normal",ifelse(ident%in%c(12,9,6,8),"LT","ST"))
seurat$status<-status

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$status
seurat<-subset(seurat,idents=c("Normal","ST"))

Idents(seurat)<-seurat$UMAP_clusters
seurat<-RenameIdents(seurat,"1"="Naive","5"="Naive","7"="Naive",
		     "2"="Memory","4"="Memory",
		     "Pb"="3","6"="Immature")

seurat$celltype<-Idents(seurat)
genes<-c("IL4R","fcn1","mzb1","igha1","CD38","irf1","ighm","ighg1","ighg2","ighv3-23","ighv3-21","ighv3-13","igkv3-11","igkv3-20")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, features =gene, split.by = "status", group.by = "celltype", pt.size = 0)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

