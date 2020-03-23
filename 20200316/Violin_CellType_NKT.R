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
seurat<-RenameIdents(seurat,"4"="CD4","5"="CD4","8"="CD4","9"="CD4",
		     "10"="CD4","12"="CD4","14"="CD4","20"="CD4","21"="CD4",
		     "2"="CD8","6"="CD8","7"="CD8","11"="CD8","13"="CD8",
		     "15"="CD8","16"="CD8","17"="CD8",
		     "1"="NK","3"="NK","19"="NK",
		     "18"="Proliferation")
seurat$celltype<-Idents(seurat)
genes<-c("JUN","FOS","JUNB","KLF6","CCL4","CCL5","CXCR4","IFRD1","IRF1","IFI6","CD38","isg15","S100A8","S100A9","STAT3","IFNg","MKI67","TYMS","FASLG")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, features =gene, split.by = "status", group.by = "celltype", pt.size = 0)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

