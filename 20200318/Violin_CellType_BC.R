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

genes<-rownames(seurat)
drop_genes<-c("IFITM1","FTL","FTH1","CD74","MALAT1")
genes<-setdiff(genes,drop_genes)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

print("Subset seurat")
Idents(seurat)<-seurat$status
seurat<-subset(seurat,idents=c("HC","ER"))
seurat$status<-factor(seurat$status,levels=c("HC","ER"))
Idents(seurat)<-seurat$UMAP_clusters
seurat<-RenameIdents(seurat,"1"="Naive","6"="Memory","7"="Naive",
		     "2"="Memory","4"="Pb","3"="Memory",
		     "5"="Pb","8"="Immature")


seurat$celltype<-Idents(seurat)
Idents(seurat)<-seurat$status
genes<-c("IL4R","fcn1","mzb1","igha1","CD38","irf1","ighm","ighg1","ighg2","ighv3-23","ighv3-21","ighv3-13",
	 "igkv3-11","igkv3-20")
genes<-str_to_upper(genes)

for(gene in genes){
	p<-VlnPlot(seurat, features =gene, split.by = "status", group.by = "celltype", pt.size = 0)+
		theme(axis.text.x=element_text(size=22),
		      plot.title = element_text(size=25),
		      legend.text=element_text(size=22),
		      axis.title.x=element_blank())
	jpeg(file.path(args$outdir,paste0("HC_ER_",gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

