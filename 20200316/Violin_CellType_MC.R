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
seurat<-RenameIdents(seurat,"9"="DC","12"="DC","13"="DC","16"="DC",
		     "2"="CD14","3"="CD14","4"="CD14","6"="CD14","7"="CD14","14"="CD14",
		     "5"="CD16","15"="CD16","17"="CD16",
		     "1"="Med","8"="Med","10"="Med","11"="Med")
seurat$celltype<-Idents(seurat)
genes<-c("JUN","junb","ifi6","tnfrsf1a","FOS","JUNB","ccl3","cxcl8","cxcr4","klf4","klf6","IL1B","CD83","irf1","ccl4","S100A9","S100A8","isg15")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, features =gene, split.by = "status", group.by = "celltype", pt.size = 0)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

