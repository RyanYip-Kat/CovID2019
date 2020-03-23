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
seurat<-RenameIdents(seurat,"1"="NK","5"="NK",
                     "10"="CD4","3"="CD4","4"="CD4",
                     "7"="CD4","9"="CD4","16"="CD4",
                     "11"="CD8","12"="CD8","13"="CD8",
                     "14"="CD8","2"="CD8","6"="CD8","8"="CD8",
                     "15"="Proliferation")

seurat$celltype<-Idents(seurat)
genes<-c("JUN","FOS","JUNB","KLF6","CCL4","CCL5","CXCR4","IFRD1","IRF1","IFI6","CD38","isg15","S100A8","S100A9","STAT3","IFNg","MKI67","TYMS","FASLG")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, features =gene, split.by = "status", group.by = "celltype", pt.size = 0)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

